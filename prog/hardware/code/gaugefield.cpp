#include "gaugefield.hpp"

#include <cmath>

#include "../../logger.hpp"
#include "../../meta/util.hpp"
#include "../device.hpp"

using namespace std;

static std::string collect_build_options(hardware::Device * device, const meta::Inputparameters& params);

const hardware::buffers::SU3 * hardware::code::Gaugefield::get_gaugefield()
{
	return &gaugefield;
}

static std::string collect_build_options(hardware::Device * device, const meta::Inputparameters& params)
{
	using namespace hardware::buffers;

	std::ostringstream options;
	options << "-D _INKERNEL_ -D NSPACE=" << params.get_nspace() << " -D NTIME=" << params.get_ntime() << " -D VOLSPACE=" << meta::get_volspace(params) << " -D VOL4D=" << meta::get_vol4d(params);

	//this is needed for hmc_ocl_su3matrix
	options << " -D SU3SIZE=" << NC*NC << " -D STAPLEMATRIXSIZE=" << NC*NC;

	if(params.get_precision() == 64) {
		options << " -D _USEDOUBLEPREC_";
		// TODO renable support for older AMD GPUs
		//if( device_double_extension.empty() ) {
		//  logger.warn() << "Warning: Undefined extension for use of double.";
		//} else {
		//  options << " -D _DEVICE_DOUBLE_EXTENSION_" << device_double_extension << "_";
		//}
		options << " -D _DEVICE_DOUBLE_EXTENSION_KHR_";
	}
	if( device->get_device_type() == CL_DEVICE_TYPE_GPU )
		options << " -D _USEGPU_";
	if(params.get_use_chem_pot_re() == true) {
		options << " -D _CP_REAL_";
		options << " -D CPR=" << params.get_chem_pot_re();
		options << " -D EXPCPR=" << exp(params.get_chem_pot_re() );
		options << " -D MEXPCPR=" << exp(-1.*params.get_chem_pot_re() );
	}
	if(params.get_use_chem_pot_im() == true) {
		options << " -D _CP_IMAG_";
		options << " -D CPI=" << params.get_chem_pot_im();
		options << " -D COSCPI=" << cos( params.get_chem_pot_im() );
		options << " -D SINCPI=" << sin( params.get_chem_pot_im() );
	}
	if(params.get_use_smearing() == true) {
		options << " -D _USE_SMEARING_";
		options << " -D RHO=" << params.get_rho();
		options << " -D RHO_ITER=" << params.get_rho_iter();
	}
	if(device->get_prefers_soa()) {
		options << " -D _USE_SOA_";
	}
	if(check_SU3_for_SOA(device)) {
		options << " -D GAUGEFIELD_STRIDE=" << get_SU3_buffer_stride(meta::get_vol4d(params) * NDIM, device);
	}
	options << " -I " << SOURCEDIR;

	if(device->get_prefers_blocked_loops()) {
		options << " -D _USE_BLOCKED_LOOPS_";
	}

	if(meta::get_use_rectangles(params) == true) {
		options <<  " -D _USE_RECT_" ;
	}
	if(params.get_use_rec12() == true) {
		options <<  " -D _USE_REC12_" ;
	}

	return options.str();
}

void hardware::code::Gaugefield::fill_kernels()
{
	basic_opencl_code = ClSourcePackage(collect_build_options(get_device(), get_parameters()))
	                    << "opencl_header.cl" << "operations_geometry.cl" << "operations_complex.cl"
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

void hardware::code::Gaugefield::clear_kernels()
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

void hardware::code::Gaugefield::plaquette_device(const hardware::buffers::SU3 * gf, const hardware::buffers::Plain<hmc_float> * plaq, const hardware::buffers::Plain<hmc_float> * tplaq, const hardware::buffers::Plain<hmc_float> * splaq)
{
	using namespace hardware::buffers;

	//query work-sizes for kernel
	size_t ls, gs;
	cl_uint num_groups;
	this->get_work_sizes(plaquette, &ls, &gs, &num_groups);

	const Plain<hmc_float> clmem_plaq_buf_glob(num_groups, get_device());
	const Plain<hmc_float> clmem_tplaq_buf_glob(num_groups, get_device());
	const Plain<hmc_float> clmem_splaq_buf_glob(num_groups, get_device());

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

	get_device()->enqueue_kernel(plaquette, gs, ls);

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
	get_device()->enqueue_kernel(plaquette_reduction, gs, ls);
}

void hardware::code::Gaugefield::rectangles_device(const hardware::buffers::SU3 * gf, const hardware::buffers::Plain<hmc_float> * rect)
{
	//query work-sizes for kernel
	size_t ls, gs;
	cl_uint num_groups;
	this->get_work_sizes(rectangles, &ls, &gs, &num_groups);

	const hardware::buffers::Plain<hmc_float> clmem_rect_buf_glob(num_groups, get_device());

	int buf_loc_size_float = sizeof(hmc_float) * ls;

	//set arguments
	// run local rectangles calculation and first part of reduction
	int clerr = clSetKernelArg(rectangles, 0, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(rectangles, 1, sizeof(cl_mem), clmem_rect_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(rectangles, 2, buf_loc_size_float, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(rectangles, gs, ls);

	// run second part of rectangles reduction

	this->get_work_sizes(rectangles_reduction, &ls, &gs, &num_groups);

	clerr = clSetKernelArg(rectangles_reduction, 0, sizeof(cl_mem), clmem_rect_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(rectangles_reduction, 1, sizeof(cl_mem), rect->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(rectangles_reduction, 2, sizeof(cl_uint), &num_groups);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);


	///@todo improve
	ls = 1;
	gs = 1;
	get_device()->enqueue_kernel(rectangles_reduction, gs, ls);

}

void hardware::code::Gaugefield::polyakov_device(const hardware::buffers::SU3 * gf, const hardware::buffers::Plain<hmc_complex> * pol)
{
	//query work-sizes for kernel
	size_t ls, gs;
	cl_uint num_groups;
	this->get_work_sizes(polyakov, &ls, &gs, &num_groups);
	int buf_loc_size_complex = sizeof(hmc_complex) * ls;

	const hardware::buffers::Plain<hmc_complex> clmem_polyakov_buf_glob(num_groups, get_device());

	// local polyakov compuation and first part of reduction
	int clerr = clSetKernelArg(polyakov, 0, sizeof(cl_mem), gf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(polyakov, 1, sizeof(cl_mem), clmem_polyakov_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(polyakov, 2, buf_loc_size_complex, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(polyakov, gs, ls);

	// second part of polyakov reduction

	this->get_work_sizes(polyakov_reduction, &ls, &gs, &num_groups);

	clerr = clSetKernelArg(polyakov_reduction, 0, sizeof(cl_mem), clmem_polyakov_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(polyakov_reduction, 1, sizeof(cl_mem), pol->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(polyakov_reduction, 2, sizeof(cl_uint), &num_groups);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	///@todo improve
	ls = 1;
	gs = 1;
	get_device()->enqueue_kernel(polyakov_reduction, gs, ls);

}

void hardware::code::Gaugefield::gaugeobservables(hmc_float * plaq_out, hmc_float * tplaq_out, hmc_float * splaq_out, hmc_complex * pol_out)
{
	gaugeobservables(get_gaugefield(), plaq_out, tplaq_out, splaq_out, pol_out);
}

void hardware::code::Gaugefield::gaugeobservables(const hardware::buffers::SU3 * gf, hmc_float * plaq_out, hmc_float * tplaq_out, hmc_float * splaq_out, hmc_complex * pol_out)
{
	const hardware::buffers::Plain<hmc_float> plaq(1, get_device());
	const hardware::buffers::Plain<hmc_float> splaq(1, get_device());
	const hardware::buffers::Plain<hmc_float> tplaq(1, get_device());
	const hardware::buffers::Plain<hmc_complex> pol(1, get_device());

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
	polyakov_device(gf, &pol);

	//read out values
	hmc_complex tmp_pol = hmc_complex_zero;
	//NOTE: this is a blocking call!
	pol.dump(&tmp_pol);

	tmp_pol.re /= static_cast<hmc_float> ( meta::get_poly_norm(get_parameters()) );
	tmp_pol.im /= static_cast<hmc_float> ( meta::get_poly_norm(get_parameters()) );

	pol_out->re = tmp_pol.re;
	pol_out->im = tmp_pol.im;
}

void hardware::code::Gaugefield::gaugeobservables_rectangles(const hardware::buffers::SU3 * gf, hmc_float * rect_out)
{
	const hardware::buffers::Plain<hmc_float> rect(1, get_device());

	//measure plaquette
	rectangles_device(gf, &rect);

	//read out values
	//NOTE: these are blocking calls!
	rect.dump(rect_out);
	//NOTE: the rectangle value has not been normalized since it is mostly used for the HMC where one needs the absolute value
}

void hardware::code::Gaugefield::stout_smear_device(const hardware::buffers::SU3 * in, const hardware::buffers::SU3 * out)
{
	//query work-sizes for kernel
	size_t ls, gs;
	cl_uint num_groups;
	this->get_work_sizes(stout_smear, &ls, &gs, &num_groups);

	int clerr = clSetKernelArg(stout_smear, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(stout_smear, 1, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(stout_smear , gs, ls);

	return;
}

void hardware::code::Gaugefield::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
	Opencl_Module::get_work_sizes(kernel, ls, gs, num_groups);
}

size_t hardware::code::Gaugefield::get_read_write_size(const std::string& in) const
{
	//Depending on the compile-options, one has different sizes...
	size_t D = meta::get_float_size(get_parameters());
	size_t R = meta::get_mat_size(get_parameters());
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

uint64_t hardware::code::Gaugefield::get_flop_size(const std::string& in) const
{
	const size_t VOL4D = meta::get_vol4d(get_parameters());
	const size_t VOLSPACE = meta::get_volspace(get_parameters());
	if (in == "polyakov") {
		//this kernel performs NTIME -1 su3matrix-multiplications, takes a complex trace and adds these real values over VOLSPACE
		return VOLSPACE * ( (get_parameters().get_ntime() - 1) * meta::get_flop_su3_su3() + meta::get_flop_su3trace()) ;
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

void hardware::code::Gaugefield::print_profiling(const std::string& filename, int number) const
{
	Opencl_Module::print_profiling(filename, number);
	Opencl_Module::print_profiling(filename, polyakov);
	Opencl_Module::print_profiling(filename, polyakov_reduction);
	Opencl_Module::print_profiling(filename, plaquette);
	Opencl_Module::print_profiling(filename, rectangles);
	Opencl_Module::print_profiling(filename, plaquette_reduction);
	Opencl_Module::print_profiling(filename, stout_smear);
	Opencl_Module::print_profiling(filename, convertGaugefieldToSOA);
	Opencl_Module::print_profiling(filename, convertGaugefieldFromSOA);
}

void hardware::code::Gaugefield::smear_gaugefield(const hardware::buffers::SU3 * gf, const std::vector<const hardware::buffers::SU3*>& gf_intermediate)
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
		hardware::buffers::SU3 gf_tmp(gf->get_elements(), get_device());
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

void hardware::code::Gaugefield::unsmear_gaugefield(const hardware::buffers::SU3 * gf)
{
	logger.debug() << "\t\trestore unsmeared gaugefield...";
	hardware::buffers::copyData(gf, &gf_unsmeared);
}

void hardware::code::Gaugefield::importGaugefield(const Matrixsu3 * const data)
{
	importGaugefield(get_gaugefield(), data);
}
void hardware::code::Gaugefield::importGaugefield(const hardware::buffers::SU3 * gaugefield, const Matrixsu3 * const data)
{
	using namespace hardware::buffers;

	logger.trace() << "Import gaugefield to get_device()";
	if(get_device()->get_prefers_soa()) {
		Plain<Matrixsu3> tmp(gaugefield->get_elements(), get_device());
		tmp.load(data);
		convertGaugefieldToSOA_device(gaugefield, &tmp);
	} else {
		gaugefield->load(data);
	}
}

void hardware::code::Gaugefield::exportGaugefield(Matrixsu3 * const dest)
{
	exportGaugefield(dest, &gaugefield);
}

void hardware::code::Gaugefield::exportGaugefield(Matrixsu3 * const dest, const hardware::buffers::SU3 * gaugefield)
{
	using namespace hardware::buffers;

	logger.trace() << "Exporting gaugefield from get_device()";
	if(get_device()->get_prefers_soa()) {
		Plain<Matrixsu3> tmp(gaugefield->get_elements(), get_device());
		convertGaugefieldFromSOA_device(&tmp, gaugefield);
		tmp.dump(dest);
	} else {
		gaugefield->dump(dest);
	}
}

void hardware::code::Gaugefield::convertGaugefieldToSOA_device(const hardware::buffers::SU3 * out, const hardware::buffers::Plain<Matrixsu3> * in)
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

	get_device()->enqueue_kernel(convertGaugefieldToSOA, gs2, ls2);
}

void hardware::code::Gaugefield::convertGaugefieldFromSOA_device(const hardware::buffers::Plain<Matrixsu3> * out, const hardware::buffers::SU3 * in)
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

	get_device()->enqueue_kernel(convertGaugefieldFromSOA, gs2, ls2);
}

hardware::code::Gaugefield::Gaugefield(const meta::Inputparameters& params, hardware::Device * device)
	: Opencl_Module(params, device), gaugefield(NDIM * meta::get_vol4d(params), device),
	  gf_unsmeared(gaugefield.get_elements(), device),
	  stout_smear(0), rectangles(0), rectangles_reduction(0)
{
	fill_kernels();
};


hardware::code::Gaugefield::~Gaugefield()
{
	clear_kernels();
}

ClSourcePackage hardware::code::Gaugefield::get_sources() const noexcept
{
	return basic_opencl_code;
}
