#include "gaugemomentum.hpp"

#include "../../logger.hpp"
#include "../../meta/util.hpp"
#include "../device.hpp"
#include <fstream>
#include <cmath>
#include "fermions.hpp"
#include "spinors.hpp"
#include "prng.hpp"

using namespace std;

static std::string collect_build_options(hardware::Device * device, const meta::Inputparameters& params);

hardware::code::Gaugemomentum::Gaugemomentum(const meta::Inputparameters& params, hardware::Device * device)
	: Opencl_Module(params, device)
{
	fill_kernels();
}

hardware::code::Gaugemomentum::~Gaugemomentum()
{
	clear_kernels();
}

static std::string collect_build_options(hardware::Device * device, const meta::Inputparameters& params)
{
	using namespace hardware::buffers;
	using namespace hardware::code;

	const size_4 mem_size = device->get_mem_lattice_size();
	const size_4 local_size = device->get_local_lattice_size();

	std::ostringstream options;
	options <<  " -D GAUGEMOMENTASIZE_GLOBAL=" << meta::get_vol4d(params) * NDIM;
	options <<  " -D GAUGEMOMENTASIZE_LOCAL=" << get_vol4d(local_size) * NDIM;
	options <<  " -D GAUGEMOMENTASIZE_MEM=" << get_vol4d(mem_size) * NDIM;
	//in case of tlsym gauge action
	if(meta::get_use_rectangles(params) == true) {
		options <<  " -D C0=" << meta::get_c0(params) << " -D C1=" << meta::get_c1(params);
	}
	if(check_Gaugemomentum_for_SOA(device)) {
		options << " -D GAUGEMOMENTA_STRIDE=" << get_Gaugemomentum_buffer_stride(get_vol4d(mem_size) * NDIM, device);
	}
	return options.str();
}

void hardware::code::Gaugemomentum::fill_kernels()
{
	basic_gaugemomentum_code = get_device()->get_fermion_code()->get_sources() << ClSourcePackage(collect_build_options(get_device(), get_parameters())) << "types_hmc.h" << "operations_gaugemomentum.cl";
	ClSourcePackage prng_code = get_device()->get_prng_code()->get_sources();

	_set_zero_gaugemomentum = createKernel("set_zero_gaugemomentum") << basic_gaugemomentum_code <<  "gaugemomentum_zero.cl";
	generate_gaussian_gaugemomenta = createKernel("generate_gaussian_gaugemomenta") << basic_gaugemomentum_code << prng_code << "gaugemomentum_gaussian.cl";
	gaugemomentum_squarenorm = createKernel("gaugemomentum_squarenorm") << basic_gaugemomentum_code << "gaugemomentum_squarenorm.cl";
	if(get_device()->get_prefers_soa()) {
		gaugemomentum_convert_to_soa = createKernel("gaugemomentum_convert_to_soa") << basic_gaugemomentum_code << "gaugemomentum_convert.cl";
		gaugemomentum_convert_from_soa = createKernel("gaugemomentum_convert_from_soa") << basic_gaugemomentum_code << "gaugemomentum_convert.cl";
	} else {
		gaugemomentum_convert_to_soa = 0;
		gaugemomentum_convert_from_soa = 0;
	}
}

void hardware::code::Gaugemomentum::clear_kernels()
{
	cl_uint clerr = CL_SUCCESS;
	///@todo some kernels are missing here
	logger.debug() << "release gaugemomentum kernels.." ;
	clerr = clReleaseKernel(generate_gaussian_gaugemomenta);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(_set_zero_gaugemomentum);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
}

void hardware::code::Gaugemomentum::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
	Opencl_Module::get_work_sizes(kernel, ls, gs, num_groups);

	// kernels that use random numbers must not exceed the size of the random state array
	if(kernel == generate_gaussian_gaugemomenta) {
		if(*gs > hardware::buffers::get_prng_buffer_size(get_device())) {
			*gs = hardware::buffers::get_prng_buffer_size(get_device());
		}
	} else if(kernel == gaugemomentum_squarenorm) {
		if(*ls != 1) { // avoid potential reduction error on AMD platform
			*ls = 128;
		}
	}

	return;
}

size_t hardware::code::Gaugemomentum::get_read_write_size(const std::string& in) const
{
	//Depending on the compile-options, one has different sizes...
	size_t D = meta::get_float_size(get_parameters());
	//this is the number of links in the system (and of gaugemomenta)
	size_t G = meta::get_vol4d(get_parameters()) * NDIM;
	//this is the same as in the function above
	//NOTE: 1 ae has NC*NC-1 = 8 real entries
	int A = meta::get_su3algebrasize();
	if (in == "generate_gaussian_gaugemomenta") {
		//this kernel writes 1 ae
		return (A) * D * G;
	}
	if (in == "set_zero_gaugemomentum;") {
		//this kernel writes 1 ae per link
		return G * D * A;
	}
	if (in == "gaugemomentum_squarenorm") {
		//this kernel reads 1 ae and writes 1 float per link
		return G * D * ( A + 1);
	}
	if (in == "gaugemomentum_convert_to_soa") {
		return 1000000000000000000000000;
	}
	if (in == "gaugemomentum_convert_from_soa") {
		return 1000000000000000000000000;
	}
	return 0;
}

uint64_t hardware::code::Gaugemomentum::get_flop_size(const std::string& in) const
{
	//this is the number of links in the system (and of gaugemomenta)
	uint64_t G = meta::get_vol4d(get_parameters()) * NDIM;
	//NOTE: 1 ae has NC*NC-1 = 8 real entries
	uint64_t A = meta::get_su3algebrasize();
	if (in == "generate_gaussian_gaugemomenta") {
		//this kernel performs 0 multiplications per site
		///@todo ? I did not count the gaussian normal pair production, which is very complicated...
		return 0;
	}
	if (in == "set_zero_gaugemomentum;") {
		//this kernel performs 0 mults
		return 0;
	}
	if (in == "gaugemomentum_squarenorm") {
		//this kernel performs 8 real mults and 8-1 real adds per ae
		return (8 + 7) * A * G;
	}
	if (in == "gaugemomentum_convert_to_soa") {
		return 1000000000000000000000000;
	}
	if (in == "gaugemomentum_convert_from_soa") {
		return 1000000000000000000000000;
	}
	return 0;
}

void hardware::code::Gaugemomentum::print_profiling(const std::string& filename, int number) const
{
	Opencl_Module::print_profiling(filename, number);
	Opencl_Module::print_profiling(filename, generate_gaussian_gaugemomenta);
	Opencl_Module::print_profiling(filename, _set_zero_gaugemomentum);
	Opencl_Module::print_profiling(filename, gaugemomentum_squarenorm);
	Opencl_Module::print_profiling(filename, gaugemomentum_convert_to_soa);
	Opencl_Module::print_profiling(filename, gaugemomentum_convert_from_soa);
}

/////////////////////////////////////////////////
// Methods on device

void hardware::code::Gaugemomentum::generate_gaussian_gaugemomenta_device(const hardware::buffers::Gaugemomentum * in, const hardware::buffers::PRNGBuffer * prng) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(generate_gaussian_gaugemomenta, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(generate_gaussian_gaugemomenta, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(generate_gaussian_gaugemomenta, 1, sizeof(cl_mem), prng->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( generate_gaussian_gaugemomenta , gs2, ls2);

	if(logger.beDebug()) {
		hardware::buffers::Plain<hmc_float> force_tmp(1, get_device());
		hmc_float resid;
		this->set_float_to_gaugemomentum_squarenorm_device(in, &force_tmp);
		force_tmp.dump(&resid);
		logger.debug() <<  "\tgaussian gaugemomenta:\t" << resid;
		if(resid != resid) {
			throw Print_Error_Message("calculation of gaussian gm gave nan! Aborting...", __FILE__, __LINE__);
		}
		if(resid == INFINITY) {
			bool writeout = false;
			if(writeout) {
				//create buffer to store ae-field
				int ae_num = in->get_elements();

				ae * ae_tmp = new ae[ae_num];

				//get buffer from device
				cout << "copy buffer to host" << endl;
				exportGaugemomentumBuffer(ae_tmp, in);

				//write out to file
				ofstream out("clmem_p_at_inf");
				if(!out) {
					cout << "Cannot open file.\n";
				}
				for(int i = 0; i < ae_num; i++) {
					out << i << "\t" << ae_tmp[i].e0 << endl;
					out << i << "\t" << ae_tmp[i].e1 << endl;
					out << i << "\t" << ae_tmp[i].e2 << endl;
					out << i << "\t" << ae_tmp[i].e3 << endl;
					out << i << "\t" << ae_tmp[i].e4 << endl;
					out << i << "\t" << ae_tmp[i].e5 << endl;
					out << i << "\t" << ae_tmp[i].e6 << endl;
					out << i << "\t" << ae_tmp[i].e7 << endl;
				}
				out.close();

				//calc sqnorm of ae_tmp
				hmc_float sqnorm = 0.;
				for(int i = 0; i < ae_num; i++) {
					sqnorm += ae_tmp[i].e0 * ae_tmp[i].e0;
					sqnorm += ae_tmp[i].e1 * ae_tmp[i].e1;
					sqnorm += ae_tmp[i].e2 * ae_tmp[i].e2;
					sqnorm += ae_tmp[i].e3 * ae_tmp[i].e3;
					sqnorm += ae_tmp[i].e4 * ae_tmp[i].e4;
					sqnorm += ae_tmp[i].e5 * ae_tmp[i].e5;
					sqnorm += ae_tmp[i].e6 * ae_tmp[i].e6;
					sqnorm += ae_tmp[i].e7 * ae_tmp[i].e7;
				}
				cout << "sqnrom: " << sqnorm << endl;
				free(ae_tmp);
			}
			throw Print_Error_Message("calculation of gaussian gm gave inf! Aborting...", __FILE__, __LINE__);
		}

	}

}

void hardware::code::Gaugemomentum::set_zero_gaugemomentum(const hardware::buffers::Gaugemomentum * buf) const
{
#ifdef CL_VERSION_1_2
	buf->clear();
#else
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(_set_zero_gaugemomentum, &ls2, &gs2, &num_groups);
	//set arguments
	//this is always applied to clmem_force
	int clerr = clSetKernelArg(_set_zero_gaugemomentum, 0, sizeof(cl_mem), buf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(_set_zero_gaugemomentum , gs2, ls2);
#endif
}


void hardware::code::Gaugemomentum::set_float_to_gaugemomentum_squarenorm_device(const hardware::buffers::Gaugemomentum * clmem_in, const hardware::buffers::Plain<hmc_float> * out) const
{
	auto spinor_code = get_device()->get_spinor_code();

	//__kernel void gaugemomentum_squarenorm(__global ae * in, __global hmc_float * out){
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(gaugemomentum_squarenorm, &ls2, &gs2, &num_groups);

	const hardware::buffers::Plain<hmc_float> clmem_global_squarenorm_buf_glob(num_groups, get_device());

	int clerr = clSetKernelArg(gaugemomentum_squarenorm, 0, sizeof(cl_mem), clmem_in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(gaugemomentum_squarenorm, 1, sizeof(cl_mem), clmem_global_squarenorm_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(gaugemomentum_squarenorm, 2, sizeof(hmc_float) * ls2, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	get_device()->enqueue_kernel(gaugemomentum_squarenorm, gs2, ls2);

	spinor_code->global_squarenorm_reduction(out, &clmem_global_squarenorm_buf_glob);
}

void hardware::code::Gaugemomentum::importGaugemomentumBuffer(const hardware::buffers::Gaugemomentum * dest, const ae * const data) const
{
	size_t const REQUIRED_BUFFER_SIZE = get_vol4d(get_device()->get_mem_lattice_size()) * NDIM;
	if(dest->get_elements() != REQUIRED_BUFFER_SIZE) {
		throw std::invalid_argument("Destination buffer is not of proper size");
	}

	cl_int clerr;
	if(dest->is_soa()) {
		hardware::buffers::Plain<ae> tmp(REQUIRED_BUFFER_SIZE, dest->get_device());
		tmp.load(data);

		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(gaugemomentum_convert_to_soa, &ls2, &gs2, &num_groups);
		clerr = clSetKernelArg(gaugemomentum_convert_to_soa, 0, sizeof(cl_mem), dest->get_cl_buffer());
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		clerr = clSetKernelArg(gaugemomentum_convert_to_soa, 1, sizeof(cl_mem), tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		get_device()->enqueue_kernel(gaugemomentum_convert_to_soa, gs2, ls2);
	} else {
		dest->load(data);
	}
}

void hardware::code::Gaugemomentum::exportGaugemomentumBuffer(ae * const dest, const hardware::buffers::Gaugemomentum * buf) const
{
	size_t const REQUIRED_BUFFER_SIZE = get_vol4d(get_device()->get_mem_lattice_size()) * NDIM;
	if(buf->get_elements() != REQUIRED_BUFFER_SIZE) {
		throw std::invalid_argument("Source buffer is not of proper size");
	}

	cl_int clerr;
	if(buf->is_soa()) {
		hardware::buffers::Plain<ae> tmp(REQUIRED_BUFFER_SIZE, buf->get_device());

		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(gaugemomentum_convert_from_soa, &ls2, &gs2, &num_groups);
		clerr = clSetKernelArg(gaugemomentum_convert_from_soa, 0, sizeof(cl_mem), tmp);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		clerr = clSetKernelArg(gaugemomentum_convert_from_soa, 1, sizeof(cl_mem), buf->get_cl_buffer());
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		get_device()->enqueue_kernel(gaugemomentum_convert_from_soa, gs2, ls2);

		tmp.dump(dest);
	} else {
		buf->dump(dest);
	}
}


ClSourcePackage hardware::code::Gaugemomentum::get_sources() const noexcept
{
	return basic_gaugemomentum_code;
}
