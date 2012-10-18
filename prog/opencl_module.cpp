#include "opencl_module.h"

#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"
#include "meta/util.hpp"

using namespace std;

static void print_profile_header(const std::string& filename, int number);
/**
 * Print the profiling information of a specific kernel to a file.
 *
 * @param filename Name of file where data is appended.
 * @param kernelName Name of specific kernel.
 * @param time_total total execution time
 * @param calls_total total number of kernel calls
 * @param read_write_size number of bytes read and written by the kernel
 * @param flop_size amount of flops performed by the kernel
 */
static void print_profiling(const std::string& filename, const std::string& kernelName, const hardware::ProfilingData& data, size_t read_write_size, uint64_t flop_size, uint64_t sites);

void Opencl_Module::init()
{
	logger.debug() << "Device is " << device->get_name();

	cl_int clerr = clGetCommandQueueInfo(get_queue(), CL_QUEUE_CONTEXT, sizeof(cl_context), &ocl_context, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetCommandQueueInfo", __FILE__, __LINE__);

	// initialize memory usage tracking
	allocated_bytes = 0;
	max_allocated_bytes = 0;
	allocated_hostptr_bytes = 0;
	logger.trace() << "Initial memory usage (" << device->get_name() << "): " << allocated_bytes << " bytes - Maximum usage: " << max_allocated_bytes << " - Host backed memory: " << allocated_hostptr_bytes << " (Assuming stored gaugefield";

	this->fill_buffers();
	this->fill_kernels();
}

void Opencl_Module::finalize()
{

	this->clear_buffers();
	this->clear_kernels();
	return;
}


cl_context Opencl_Module::get_context()
{
	return ocl_context;
}


const hardware::buffers::SU3 * Opencl_Module::get_gaugefield()
{
	return &gaugefield;
}

const meta::Inputparameters& Opencl_Module::get_parameters() const noexcept
{
	return parameters;
}

hardware::Device * Opencl_Module::get_device() const noexcept
{
	return device;
}

void Opencl_Module::fill_collect_options(stringstream* collect_options)
{
	using namespace hardware::buffers;

	*collect_options << "-D_INKERNEL_ -DNSPACE=" << get_parameters().get_nspace() << " -DNTIME=" << get_parameters().get_ntime() << " -DVOLSPACE=" << meta::get_volspace(get_parameters()) << " -DVOL4D=" << meta::get_vol4d(get_parameters());

	//this is needed for hmc_ocl_su3matrix
	*collect_options << " -DSU3SIZE=" << NC*NC << " -DSTAPLEMATRIXSIZE=" << NC*NC;

	if(get_parameters().get_precision() == 64) {
		*collect_options << " -D_USEDOUBLEPREC_";
		// TODO renable support for older AMD GPUs
		//if( device_double_extension.empty() ) {
		//  logger.warn() << "Warning: Undefined extension for use of double.";
		//} else {
		//  *collect_options << " -D_DEVICE_DOUBLE_EXTENSION_" << device_double_extension << "_";
		//}
		*collect_options << " -D_DEVICE_DOUBLE_EXTENSION_KHR_";
	}
	if( device->get_device_type() == CL_DEVICE_TYPE_GPU )
		*collect_options << " -D_USEGPU_";
	if(get_parameters().get_use_chem_pot_re() == true) {
		*collect_options << " -D_CP_REAL_";
		*collect_options << " -DCPR=" << get_parameters().get_chem_pot_re();
		*collect_options << " -DEXPCPR=" << exp(get_parameters().get_chem_pot_re() );
		*collect_options << " -DMEXPCPR=" << exp(-1.*get_parameters().get_chem_pot_re() );
	}
	if(get_parameters().get_use_chem_pot_im() == true) {
		*collect_options << " -D_CP_IMAG_";
		*collect_options << " -DCPI=" << get_parameters().get_chem_pot_im();
		*collect_options << " -DCOSCPI=" << cos( get_parameters().get_chem_pot_im() );
		*collect_options << " -DSINCPI=" << sin( get_parameters().get_chem_pot_im() );
	}
	if(get_parameters().get_use_smearing() == true) {
		*collect_options << " -D_USE_SMEARING_";
		*collect_options << " -DRHO=" << get_parameters().get_rho();
		*collect_options << " -DRHO_ITER=" << get_parameters().get_rho_iter();
	}
	if(device->get_prefers_soa()) {
		*collect_options << " -D_USE_SOA_";
	}
	if(check_SU3_for_SOA(get_device())) {
		*collect_options << " -DGAUGEFIELD_STRIDE=" << get_SU3_buffer_stride(meta::get_vol4d(get_parameters()) * NDIM, device);
	}
	*collect_options << " -I" << SOURCEDIR;

	if(device->get_prefers_blocked_loops()) {
		*collect_options << " -D_USE_BLOCKED_LOOPS_";
	}

	if(meta::get_use_rectangles(get_parameters()) == true) {
		*collect_options <<  " -D_USE_RECT_" ;
	}
	if(get_parameters().get_use_rec12() == true) {
		*collect_options <<  " -D_USE_REC12_" ;
	}

	return;
}

cl_mem Opencl_Module::createBuffer(cl_mem_flags flags, size_t size)
{
	return createBuffer(flags, size, 0);
}

void Opencl_Module::markMemReleased(bool host, size_t size)
{
	if(host) {
		allocated_hostptr_bytes -= size;
	} else {
		allocated_bytes -= size;
	}
	logger.trace() << "Memory usage (" << device->get_name() << "): " << allocated_bytes << " bytes - Maximum usage: " << max_allocated_bytes << " - Host backed memory: " << allocated_hostptr_bytes;
}

struct MemObjectReleaseInfo {
	size_t bytes;
	bool host;
	Opencl_Module * module;

	MemObjectReleaseInfo(size_t bytes, bool host, Opencl_Module * module)
		: bytes(bytes), host(host), module(module) { };
};

void memObjectReleased(cl_mem, void * user_data)
{
	MemObjectReleaseInfo * release_info = static_cast<MemObjectReleaseInfo *>(user_data);
	release_info->module->markMemReleased(release_info->host, release_info->bytes);
	delete release_info;
}

cl_mem Opencl_Module::createBuffer(cl_mem_flags flags, size_t size, void * host_pointer)
{
	logger.trace() << "Allocating " << size << " bytes.";

	// create buffer object
	cl_int clerr;
	cl_mem tmp = clCreateBuffer(ocl_context, flags, size, host_pointer, &clerr);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clCreateBuffer", __FILE__, __LINE__);

	// take care of memory usage bookkeeping
	bool host = (flags & CL_MEM_ALLOC_HOST_PTR) || (flags & CL_MEM_USE_HOST_PTR);
	if(host) {
		allocated_hostptr_bytes += size;
	} else {
		allocated_bytes += size;
	}
	MemObjectReleaseInfo * releaseInfo = new MemObjectReleaseInfo(size, host, this);
	clerr = clSetMemObjectDestructorCallback(tmp, memObjectReleased, releaseInfo);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clCreateMemObjectDestructorCallback", __FILE__, __LINE__);
	if(allocated_bytes >= max_allocated_bytes) {
		max_allocated_bytes = allocated_bytes;
	}

	logger.trace() << "Memory usage (" << device->get_name() << "): " << allocated_bytes << " bytes - Maximum usage: " << max_allocated_bytes << " - Host backed memory: " << allocated_hostptr_bytes;

	return tmp;
}


cl_mem Opencl_Module::create_rw_buffer(size_t size)
{
	return createBuffer(CL_MEM_READ_WRITE, size);
}

cl_mem Opencl_Module::create_wo_buffer(size_t size)
{
	return createBuffer(CL_MEM_WRITE_ONLY, size);
}

cl_mem Opencl_Module::create_ro_buffer(size_t size)
{
	return createBuffer(CL_MEM_READ_ONLY, size);
}

cl_mem Opencl_Module::create_uhp_buffer(size_t size, void *host_pointer)
{
	return createBuffer(CL_MEM_USE_HOST_PTR, size, host_pointer);
}

cl_mem Opencl_Module::create_ahp_buffer(size_t size)
{
	return createBuffer(CL_MEM_ALLOC_HOST_PTR, size);
}

cl_mem Opencl_Module::create_chp_buffer(size_t size, void *host_pointer)
{
	return createBuffer(CL_MEM_COPY_HOST_PTR, size, host_pointer);
}

void Opencl_Module::fill_buffers()
{
	// FIXME to be removed
}

void Opencl_Module::fill_kernels()
{
	basic_opencl_code = ClSourcePackage() << "opencl_header.cl" << "operations_geometry.cl" << "operations_complex.cl"
	                    << "operations_matrix_su3.cl" << "operations_matrix.cl" << "operations_gaugefield.cl";

	logger.debug() << "Create gaugeobservables kernels...";
	plaquette = createKernel("plaquette") << basic_opencl_code << "gaugeobservables_plaquette.cl";
	plaquette_reduction = createKernel("plaquette_reduction") << basic_opencl_code << "gaugeobservables_plaquette.cl";
	if(meta::get_use_rectangles(get_parameters()) == true) {
		rectangles = createKernel("rectangles") << basic_opencl_code << "gaugeobservables_rectangles.cl";
		rectangles_reduction = createKernel("rectangles_reduction") << basic_opencl_code << "gaugeobservables_rectangles.cl";
	}
	polyakov = createKernel("polyakov") << basic_opencl_code << "gaugeobservables_polyakov.cl";
	polyakov_reduction = createKernel("polyakov_reduction") << basic_opencl_code << "gaugeobservables_polyakov.cl";
	if(get_parameters().get_use_smearing() == true) {
		stout_smear = createKernel("stout_smear") << basic_opencl_code << "operations_gaugemomentum.cl" << "stout_smear.cl";
	}
	convertGaugefieldToSOA = createKernel("convertGaugefieldToSOA") << basic_opencl_code << "gaugefield_convert.cl";
	convertGaugefieldFromSOA = createKernel("convertGaugefieldFromSOA") << basic_opencl_code << "gaugefield_convert.cl";
}

void Opencl_Module::clear_buffers()
{
	logger.info() << "Maximum memory used (" << device->get_name() << "): " << max_allocated_bytes << " bytes";
}

void Opencl_Module::clear_kernels()
{
	logger.trace() << "Clearing kernels";

	cl_int clerr = CL_SUCCESS;

	clerr = clReleaseKernel(plaquette);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(polyakov);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(plaquette_reduction);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(polyakov_reduction);
	if(meta::get_use_rectangles(get_parameters()) == true) {
		clerr = clReleaseKernel(rectangles);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(rectangles_reduction);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	if(get_parameters().get_use_smearing() == true) {
		clerr = clReleaseKernel(stout_smear);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	clerr = clReleaseKernel(convertGaugefieldToSOA);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(convertGaugefieldFromSOA);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
}

void Opencl_Module::copy_buffer_on_device(cl_mem in, cl_mem out, size_t size)
{
	(*this->get_copy_on()).reset();

	int clerr = clEnqueueCopyBuffer(get_queue(), in, out, 0, 0, size , 0, 0, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueCopyBuffer", __FILE__, __LINE__);

	(*this->get_copy_on()).add();
}

void Opencl_Module::copy_buffer_to_device(void * source, cl_mem dest, size_t size)
{
	(*this->get_copy_to()).reset();

	int clerr = clEnqueueWriteBuffer(get_queue(), dest, CL_TRUE, 0, size, source, 0, 0, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueWriteBuffer", __FILE__, __LINE__);

	(*this->get_copy_to()).add();
}

void Opencl_Module::get_buffer_from_device(cl_mem source, void * dest, size_t size)
{
	(*this->get_copy_to()).reset();
	cl_int clerr = clEnqueueReadBuffer(get_queue(), source, CL_TRUE, 0, size, dest, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueReadBuffer", __FILE__, __LINE__);

	(*this->get_copy_to()).add();
}

void Opencl_Module::plaquette_device(const hardware::buffers::SU3 * gf, const hardware::buffers::Plain<hmc_float> * plaq, const hardware::buffers::Plain<hmc_float> * tplaq, const hardware::buffers::Plain<hmc_float> * splaq)
{
	using namespace hardware::buffers;

	//query work-sizes for kernel
	size_t ls, gs;
	cl_uint num_groups;
	this->get_work_sizes(plaquette, &ls, &gs, &num_groups);

	const Plain<hmc_float> clmem_plaq_buf_glob(num_groups, device);
	const Plain<hmc_float> clmem_tplaq_buf_glob(num_groups, device);
	const Plain<hmc_float> clmem_splaq_buf_glob(num_groups, device);

	int buf_loc_size_float = sizeof(hmc_float) * ls;

	//set arguments
	// run local plaquette calculation and first part of reduction
	int clerr = clSetKernelArg(plaquette, 0, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(plaquette, 1, sizeof(cl_mem), clmem_plaq_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(plaquette, 2, sizeof(cl_mem), clmem_tplaq_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(plaquette, 3, sizeof(cl_mem), clmem_splaq_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(plaquette, 4, buf_loc_size_float, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(plaquette, 5, buf_loc_size_float, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(plaquette, 6, buf_loc_size_float, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	device->enqueue_kernel(plaquette, gs, ls);

	// run second part of plaquette reduction

	this->get_work_sizes(plaquette_reduction, &ls, &gs, &num_groups);

	clerr = clSetKernelArg(plaquette_reduction, 0, sizeof(cl_mem), clmem_plaq_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(plaquette_reduction, 1, sizeof(cl_mem), clmem_tplaq_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(plaquette_reduction, 2, sizeof(cl_mem), clmem_splaq_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(plaquette_reduction, 3, sizeof(cl_mem), plaq->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(plaquette_reduction, 4, sizeof(cl_mem), tplaq->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(plaquette_reduction, 5, sizeof(cl_mem), splaq->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(plaquette_reduction, 6, sizeof(cl_uint), &num_groups);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);


	///@todo improve
	ls = 1;
	gs = 1;
	device->enqueue_kernel(plaquette_reduction, gs, ls);
}

void Opencl_Module::rectangles_device(const hardware::buffers::SU3 * gf)
{
	//query work-sizes for kernel
	size_t ls, gs;
	cl_uint num_groups;
	this->get_work_sizes(rectangles, &ls, &gs, &num_groups);

	const hardware::buffers::Plain<hmc_float> clmem_rect_buf_glob(num_groups, device);

	int buf_loc_size_float = sizeof(hmc_float) * ls;

	//set arguments
	// run local rectangles calculation and first part of reduction
	int clerr = clSetKernelArg(rectangles, 0, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(rectangles, 1, sizeof(cl_mem), clmem_rect_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(rectangles, 2, buf_loc_size_float, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	device->enqueue_kernel(rectangles, gs, ls);

	// run second part of rectangles reduction

	this->get_work_sizes(rectangles_reduction, &ls, &gs, &num_groups);

	clerr = clSetKernelArg(rectangles_reduction, 0, sizeof(cl_mem), clmem_rect_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(rectangles_reduction, 1, sizeof(cl_mem), clmem_rect);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(rectangles_reduction, 2, sizeof(cl_uint), &num_groups);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);


	///@todo improve
	ls = 1;
	gs = 1;
	device->enqueue_kernel(rectangles_reduction, gs, ls);

}

void Opencl_Module::polyakov_device(const hardware::buffers::SU3 * gf)
{
	//query work-sizes for kernel
	size_t ls, gs;
	cl_uint num_groups;
	this->get_work_sizes(polyakov, &ls, &gs, &num_groups);
	int buf_loc_size_complex = sizeof(hmc_complex) * ls;

	const hardware::buffers::Plain<hmc_complex> clmem_polyakov_buf_glob(num_groups, device);

	// local polyakov compuation and first part of reduction
	int clerr = clSetKernelArg(polyakov, 0, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(polyakov, 1, sizeof(cl_mem), clmem_polyakov_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(polyakov, 2, buf_loc_size_complex, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	device->enqueue_kernel(polyakov, gs, ls);

	// second part of polyakov reduction

	this->get_work_sizes(polyakov_reduction, &ls, &gs, &num_groups);

	clerr = clSetKernelArg(polyakov_reduction, 0, sizeof(cl_mem), clmem_polyakov_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(polyakov_reduction, 1, sizeof(cl_mem), clmem_polyakov);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(polyakov_reduction, 2, sizeof(cl_uint), &num_groups);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	///@todo improve
	ls = 1;
	gs = 1;
	device->enqueue_kernel(polyakov_reduction, gs, ls);

}

void Opencl_Module::gaugeobservables(hmc_float * plaq_out, hmc_float * tplaq_out, hmc_float * splaq_out, hmc_complex * pol_out)
{
	gaugeobservables(get_gaugefield(), plaq_out, tplaq_out, splaq_out, pol_out);
}

void Opencl_Module::gaugeobservables(const hardware::buffers::SU3 * gf, hmc_float * plaq_out, hmc_float * tplaq_out, hmc_float * splaq_out, hmc_complex * pol_out)
{
	const hardware::buffers::Plain<hmc_float> plaq(1, get_device());
	const hardware::buffers::Plain<hmc_float> splaq(1, get_device());
	const hardware::buffers::Plain<hmc_float> tplaq(1, get_device());

	//measure plaquette
	plaquette_device(gf, &plaq, &tplaq, &splaq);

	//read out values
	hmc_float tmp_plaq = 0.;
	hmc_float tmp_splaq = 0.;
	hmc_float tmp_tplaq = 0.;
	//NOTE: these are blocking calls!
	plaq.dump(&tmp_plaq);
	splaq.dump(&tmp_splaq);
	tplaq.dump(&tmp_tplaq);

	tmp_tplaq /= static_cast<hmc_float>(meta::get_tplaq_norm(get_parameters()));
	tmp_splaq /= static_cast<hmc_float>(meta::get_splaq_norm(get_parameters()));
	tmp_plaq  /= static_cast<hmc_float>(meta::get_plaq_norm(get_parameters()));

	(*plaq_out) = tmp_plaq;
	(*splaq_out) = tmp_splaq;
	(*tplaq_out) = tmp_tplaq;

	//measure polyakovloop
	polyakov_device(gf);

	//read out values
	hmc_complex pol = hmc_complex_zero;
	//NOTE: this is a blocking call!
	clmem_polyakov.dump(&pol);

	pol.re /= static_cast<hmc_float> ( meta::get_poly_norm(get_parameters()) );
	pol.im /= static_cast<hmc_float> ( meta::get_poly_norm(get_parameters()) );

	pol_out->re = pol.re;
	pol_out->im = pol.im;
}

void Opencl_Module::gaugeobservables_rectangles(const hardware::buffers::SU3 * gf, hmc_float * rect_out)
{
	//measure plaquette
	rectangles_device(gf);

	//read out values
	//NOTE: these are blocking calls!
	clmem_rect.dump(rect_out);
	//NOTE: the rectangle value has not been normalized since it is mostly used for the HMC where one needs the absolute value
}

TmpClKernel Opencl_Module::createKernel(const char * const kernel_name, const char * const build_opts)
{
	stringstream collect_options;
	if(get_parameters().is_ocl_compiler_opt_disabled()) {
		collect_options << "-cl-opt-disable ";
	}
	if(build_opts) {
		collect_options << build_opts << ' ';
	}
	this->fill_collect_options(&collect_options);
	return device->create_kernel(kernel_name, collect_options.str());
}

void Opencl_Module::stout_smear_device(const hardware::buffers::SU3 * in, const hardware::buffers::SU3 * out)
{
	//query work-sizes for kernel
	size_t ls, gs;
	cl_uint num_groups;
	this->get_work_sizes(stout_smear, &ls, &gs, &num_groups);

	int clerr = clSetKernelArg(stout_smear, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(stout_smear, 1, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	device->enqueue_kernel(stout_smear , gs, ls);

	return;
}



void Opencl_Module::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
	//Query kernel name
	string kernelname = get_kernel_name(kernel);

	size_t local_work_size = device->get_preferred_local_thread_num();
	size_t global_work_size = device->get_preferred_global_thread_num();

	const cl_uint num_groups_tmp = (global_work_size + local_work_size - 1) / local_work_size;
	global_work_size = local_work_size * num_groups_tmp;

	//write out values
	*ls = local_work_size;
	*gs = global_work_size;
	*num_groups = num_groups_tmp;

	return;
}

string Opencl_Module::get_kernel_name(const cl_kernel kernel) const
{
	int clerr;
	size_t bytesInKernelName;
	clerr = clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, 0, NULL, &bytesInKernelName);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetKernelInfo", __FILE__, __LINE__);
	char * kernelName = new char[bytesInKernelName]; // additional space for terminating 0 byte
	clerr = clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, bytesInKernelName, kernelName, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clGetKernelInfo", __FILE__, __LINE__);

	string kernel_name = kernelName;
	delete [] kernelName;

	return kernel_name;
}

usetimer * Opencl_Module::get_copy_on()
{
	return &copy_on;
}

usetimer * Opencl_Module::get_copy_to()
{
	return &copy_to;
}

size_t Opencl_Module::get_read_write_size(const std::string& in) const
{
	//Depending on the compile-options, one has different sizes...
	size_t D = meta::get_float_size(parameters);
	size_t R = meta::get_mat_size(parameters);
	//factor for complex numbers
	int C = 2;
	const size_t VOL4D = meta::get_vol4d(get_parameters());
	if (in == "polyakov") {
		//this kernel reads NTIME*VOLSPACE=VOL4D su3matrices and writes NUM_GROUPS complex numbers
		//query work-sizes for kernel to get num_groups
		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(polyakov, &ls2, &gs2, &num_groups);
		return VOL4D * D * R + num_groups * C * D;
	}
	if (in == "polyakov_reduction") {
		//this kernel reads NUM_GROUPS complex numbers and writes 1 complex number
		//query work-sizes for kernel to get num_groups
		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(polyakov_reduction, &ls2, &gs2, &num_groups);
		return (num_groups + 1 ) * C * D;
	}
	if (in == "plaquette") {
		//this kernel reads in VOL4D * ND * (ND-1) su3matrices and writes 3*num_groups real numbers
		//query work-sizes for kernel to get num_groups
		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(plaquette, &ls2, &gs2, &num_groups);
		return C * 4 * 3 * VOL4D * D * R + 3 * D * num_groups;
	}
	if (in == "plaquette_reduction") {
		//this kernel reads 3*NUM_GROUPS real numbers and writes 3 real numbers
		//query work-sizes for kernel to get num_groups
		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(polyakov_reduction, &ls2, &gs2, &num_groups);
		return (num_groups + 1 ) * 3 * C * D;
	}
	if (in == "rectangles") {
		return 1000000000000000000000;
	}
	if (in == "rectangles_reduction") {
		return 1000000000000000000000;
	}
	if (in == "stout_smear") {
		//this kernel reads in a complete gaugefield + a staple on each site and writes out a complete gaugefield
		return VOL4D * NDIM * D * R * (6 * (NDIM - 1) + 1 + 1 );
	}
	if(in == "convertGaugefieldToSOA") {
		return 2 * meta::get_vol4d(get_parameters()) * NDIM * R * C * D;
	}
	if(in == "convertGaugefieldFromSOA") {
		return 2 * meta::get_vol4d(get_parameters()) * NDIM * R * C * D;
	}
	return 0;
}

uint64_t Opencl_Module::get_flop_size(const std::string& in) const
{
	const size_t VOL4D = meta::get_vol4d(get_parameters());
	const size_t VOLSPACE = meta::get_volspace(get_parameters());
	if (in == "polyakov") {
		//this kernel performs NTIME -1 su3matrix-multiplications, takes a complex trace and adds these real values over VOLSPACE
		return VOLSPACE * ( (parameters.get_ntime() - 1) * meta::get_flop_su3_su3() + meta::get_flop_su3trace()) ;
	}
	if (in == "polyakov_reduction") {
		return 1000000000000000000000;
	}
	if (in == "plaquette") {
		//this kernel performs 3 su3matrix-mutliplications, a real su3 trace and sums over VOL4D and mu and nu (nu<mu)
		return VOL4D * NDIM * (NDIM - 1) * ( 3 + NC);
	}
	if (in == "plaquette_reduction") {
		return 1000000000000000000000;
	}
	if (in == "rectangles") {
		return 1000000000000000000000;
	}
	if (in == "rectangles_reduction") {
		return 1000000000000000000000;
	}
	if (in == "stout_smear") {
		return 1000000000000000000000;
	}
	return 0;
}

static void print_profiling(const std::string& filename, const std::string& kernelName, const hardware::ProfilingData& data, size_t read_write_size, uint64_t flop_size, uint64_t sites)
{
	hmc_float bandwidth = 0.;
	hmc_float flops = 0.;
	uint64_t avg_time = 0.;
	uint64_t avg_time_site = 0.;
	//check if kernel has been called at all
	if(data.get_num_values()) {
		avg_time = (uint64_t) ( ( (float) data.get_total_time() ) / ((float) data.get_num_values() ) );
		avg_time_site = (uint64_t) ( ( (float) data.get_total_time() ) / ((float) (data.get_num_values() * sites)) );
		//Bandwidth in GB/s: 1e-3 = 1e6 (museconds) * 1e-9 (GByte)
		bandwidth = (hmc_float) read_write_size / (hmc_float) data.get_total_time() * (hmc_float) data.get_num_values() * 1e-3;
		flops = (hmc_float) flop_size / (hmc_float) data.get_total_time() * (hmc_float) data.get_num_values() * 1e-3;
	}
	float mega = 1024 * 1024;
	//write to stream
	fstream out;
	out.open(filename.c_str(), std::ios::out | std::ios::app);
	if(!out.is_open()) File_Exception(filename.c_str());
	//CP: this is set manually to fit the longest kernel name
	out.width(32);
	out.precision(15);
	//to look like that
	/*
	logger.trace() << "*******************************************************************";
	logger.trace() << "Fermion\t"<< setfill(' ') << setw(16)<< "BW[GB/s]\t" << setfill(' ') << setw(18) << "Re/Wr[MByte]\t" << setfill(' ') << setw(6)  << "Calls\t" << setfill(' ') << setw(10)  << "Time[mus]";
	*/
	out << kernelName << "\t" << data.get_total_time() << "\t" << data.get_num_values() << "\t" << avg_time << "\t" << avg_time_site << "\t" << bandwidth << "\t" << flops << "\t" << (float) read_write_size / mega << "\t" << flop_size << std::endl;
	out.close();
	return;
}

static void print_profile_header(const std::string& filename, int number)
{
	//write to stream
	fstream out;
	out.open(filename.c_str(), std::ios::out | std::ios::app);
	if(!out.is_open()) File_Exception(filename.c_str());
	//CP: this is set manually to fit the longest kernel name
	out.width(32);
	out.precision(15);
	out << "#device " << number << "\tTime [mus]\tCalls\tAvg Time [mus]\tAvg Time/Site [mus]\tBW [GB/s]\tFLOPS [GFLOP/s]\tRe/Wr [MB]\tFLOP" << std::endl;
	return;
}

void Opencl_Module::print_profiling(const std::string& filename, const cl_kernel& kernel) const
{
	// only print info if kernel has been initialized
	if(kernel) {
		const std::string kernel_name = get_kernel_name(kernel);
		const hardware::ProfilingData data = device->get_profiling_data(kernel);
		::print_profiling(filename, kernel_name, data, this->get_read_write_size(kernel_name), this->get_flop_size(kernel_name), meta::get_vol4d(get_parameters()));
	}
}

void Opencl_Module::print_profiling(const std::string& filename, int number)
{
	logger.trace() << "Printing Profiling-information to file \"" << filename << "\"";
	print_profile_header(filename, number);
	print_profiling(filename, polyakov);
	print_profiling(filename, polyakov_reduction);
	print_profiling(filename, plaquette);
	print_profiling(filename, rectangles);
	print_profiling(filename, plaquette_reduction);
	print_profiling(filename, stout_smear);
	print_profiling(filename, convertGaugefieldToSOA);
	print_profiling(filename, convertGaugefieldFromSOA);
}

void Opencl_Module::print_copy_times(uint64_t totaltime)
{
	//copy1 ^= copy_to_from_dev_time
	//copy2 ^= copy_on_dev_time

	uint64_t copy2_time = (this->copy_on).getTime();
	uint64_t copy1_time = (this->copy_to).getTime();

	int copy1_steps =  (this->copy_to).getNumMeas();
	int copy2_steps =  (this->copy_on).getNumMeas();

	uint64_t copy1_avgtime = divide(copy1_time, copy1_steps);
	uint64_t copy2_avgtime = divide(copy2_time, copy2_steps);

	logger.info() << "## *******************************************************************";
	logger.info() << "## Copy-Times\t" << setfill(' ') << setw(12) << "total" << '\t' << setw(12) << "avg" << '\t' << setw(5) << "perc";
	logger.info() << "## CpyTo:\t" << setfill(' ') << setw(12) << copy1_time << '\t' << setw(12) << copy1_avgtime << '\t' << fixed << setw(5) << setprecision(1) << percent(copy1_time, totaltime);
	logger.info() << "## CpyOn:\t" << setfill(' ') << setw(12) << copy2_time << '\t' << setw(12) << copy2_avgtime << '\t' << fixed << setw(5) << setprecision(1) << percent(copy2_time, totaltime);
	logger.info() << "## *******************************************************************";

	logger.info() << "## No output of times to file implemented yet...";
	/** @todo output to file is not implemented */
	//See older files for example code
	return;
}

void Opencl_Module::smear_gaugefield(const hardware::buffers::SU3 * gf, const std::vector<const hardware::buffers::SU3*>& gf_intermediate)
{
	logger.debug() << "\t\tsave unsmeared gaugefield...";
	// TODO what if called before
	// FIXME memory leak if not unsmeared
	hardware::buffers::copyData(&gf_unsmeared, gf);
	logger.debug() << "\t\tperform " << get_parameters().get_rho_iter() << " steps of stout-smearing to the gaugefield...";
	if(gf_intermediate.size() == get_parameters().get_rho_iter()) {
		//the first step is applied to the original gf
		stout_smear_device(gf, gf_intermediate[0]);
		//perform rho_iter -2 intermediate steps
		for(int i = 1; i < get_parameters().get_rho_iter() - 1; i++) {
			stout_smear_device(gf_intermediate[i - 1], gf_intermediate[i]);
		}
		//the last step results in the smeared gf
		stout_smear_device(gf_intermediate[get_parameters().get_rho_iter() - 1 ], gf);
	} else {
		//one needs a temporary gf to apply the smearing to
		hardware::buffers::SU3 gf_tmp(gf->get_elements(), device);
		for(int i = 0; i < get_parameters().get_rho_iter(); i++) {
			if(i / 2) stout_smear_device(gf, &gf_tmp);
			else stout_smear_device(&gf_tmp, gf);
		}
		//if rho_iter is odd one has to copy ones more
		if(get_parameters().get_rho_iter() / 2 == 1) {
			hardware::buffers::copyData(gf, &gf_tmp);
		}
	}
}

void Opencl_Module::unsmear_gaugefield(const hardware::buffers::SU3 * gf)
{
	logger.debug() << "\t\trestore unsmeared gaugefield...";
	hardware::buffers::copyData(gf, &gf_unsmeared);
}

void Opencl_Module::importGaugefield(const Matrixsu3 * const data)
{
	importGaugefield(get_gaugefield(), data);
}
void Opencl_Module::importGaugefield(const hardware::buffers::SU3 * gaugefield, const Matrixsu3 * const data)
{
	using namespace hardware::buffers;

	logger.trace() << "Import gaugefield to device";
	if(device->get_prefers_soa()) {
		Plain<Matrixsu3> tmp(gaugefield->get_elements(), device);
		tmp.load(data);
		convertGaugefieldToSOA_device(gaugefield, &tmp);
	} else {
		gaugefield->load(data);
	}
}

void Opencl_Module::exportGaugefield(Matrixsu3 * const dest)
{
	using namespace hardware::buffers;

	logger.trace() << "Exporting gaugefield from device";
	if(device->get_prefers_soa()) {
		Plain<Matrixsu3> tmp(gaugefield.get_elements(), device);
		convertGaugefieldFromSOA_device(&tmp, &gaugefield);
		tmp.dump(dest);
	} else {
		gaugefield.dump(dest);
	}
}

void Opencl_Module::convertGaugefieldToSOA_device(const hardware::buffers::SU3 * out, const hardware::buffers::Plain<Matrixsu3> * in)
{
	if(!out->is_soa()) {
		throw std::invalid_argument("Destination buffer must be a SOA buffer");
	}

	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(convertGaugefieldToSOA, &ls2, &gs2, &num_groups);

	//set arguments
	int clerr = clSetKernelArg(convertGaugefieldToSOA, 0, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(convertGaugefieldToSOA, 1, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	device->enqueue_kernel(convertGaugefieldToSOA, gs2, ls2);
}

void Opencl_Module::convertGaugefieldFromSOA_device(const hardware::buffers::Plain<Matrixsu3> * out, const hardware::buffers::SU3 * in)
{
	if(!in->is_soa()) {
		throw std::invalid_argument("Source buffer must be a SOA buffer");
	}

	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(convertGaugefieldFromSOA, &ls2, &gs2, &num_groups);

	//set arguments
	int clerr = clSetKernelArg(convertGaugefieldFromSOA, 0, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(convertGaugefieldFromSOA, 1, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	device->enqueue_kernel(convertGaugefieldFromSOA, gs2, ls2);
}

cl_command_queue Opencl_Module::get_queue() const noexcept
{
	return device->get_queue();
}
