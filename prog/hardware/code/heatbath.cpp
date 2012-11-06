#include "heatbath.hpp"

#include "../../logger.hpp"
#include "../../meta/util.hpp"
#include "../device.hpp"

using namespace std;

static std::string collect_build_options(hardware::Device * device, const meta::Inputparameters& params);

static std::string collect_build_options(hardware::Device *, const meta::Inputparameters& params)
{
	std::ostringstream options;

	options <<  "-D BETA=" << params.get_beta();
	if(params.get_use_aniso() == true) {
		options << " -D _ANISO_";
		options <<  " -D XI_0=" << meta::get_xi_0(params);
	}

	return options.str();
}

void hardware::code::Heatbath::fill_kernels()
{
	ClSourcePackage sources = get_device()->get_gaugefield_code()->get_sources()
	                          << get_device()->get_prng_code()->get_sources()
	                          << ClSourcePackage(collect_build_options(get_device(), get_parameters()));

	logger.debug() << "Create heatbath kernels...";
	heatbath_even = createKernel("heatbath_even") << sources << "operations_heatbath.cl" << "heatbath_even.cl";
	heatbath_odd = createKernel("heatbath_odd") << sources << "operations_heatbath.cl" << "heatbath_odd.cl";

	logger.debug() << "Create overrelax kernels...";
	overrelax_even = createKernel("overrelax_even") << sources << "operations_heatbath.cl" << "overrelax_even.cl";
	overrelax_odd = createKernel("overrelax_odd") << sources << "operations_heatbath.cl" << "overrelax_odd.cl";
}

void hardware::code::Heatbath::clear_kernels()
{
	cl_int clerr = CL_SUCCESS;

	clerr = clReleaseKernel(heatbath_even);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);

	clerr = clReleaseKernel(heatbath_odd);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);

	clerr = clReleaseKernel(overrelax_even);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);

	clerr = clReleaseKernel(overrelax_odd);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
}

void hardware::code::Heatbath::run_heatbath(const hardware::buffers::SU3 * gaugefield, const hardware::buffers::PRNGBuffer * prng) const
{
	cl_int clerr = CL_SUCCESS;

	size_t global_work_size, ls;
	cl_uint num_groups;
	this->get_work_sizes(heatbath_even, &ls, &global_work_size, &num_groups);

	clerr = clSetKernelArg(heatbath_even, 0, sizeof(cl_mem), gaugefield->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(heatbath_even, 2, sizeof(cl_mem), prng->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	for(cl_int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(heatbath_even, 1, sizeof(cl_int), &i);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		get_device()->enqueue_kernel(heatbath_even, global_work_size, ls);
	}

	this->get_work_sizes(heatbath_odd, &ls, &global_work_size, &num_groups);

	clerr = clSetKernelArg(heatbath_odd, 0, sizeof(cl_mem), gaugefield->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(heatbath_odd, 2, sizeof(cl_mem), prng->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	for(cl_int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(heatbath_odd, 1, sizeof(cl_int), &i);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		get_device()->enqueue_kernel(heatbath_odd, global_work_size, ls);
	}
}

void hardware::code::Heatbath::run_overrelax(const hardware::buffers::SU3 * gaugefield, const hardware::buffers::PRNGBuffer * prng) const
{
	cl_int clerr = CL_SUCCESS;

	size_t global_work_size, ls;
	cl_uint num_groups;
	this->get_work_sizes(overrelax_even, &ls, &global_work_size, &num_groups);

	clerr = clSetKernelArg(overrelax_even, 0, sizeof(cl_mem), gaugefield->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(overrelax_even, 2, sizeof(cl_mem), prng->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	for(cl_int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(overrelax_even, 1, sizeof(cl_int), &i);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		get_device()->enqueue_kernel(overrelax_even, global_work_size, ls);
	}

	this->get_work_sizes(overrelax_odd, &ls, &global_work_size, &num_groups);

	clerr = clSetKernelArg(overrelax_odd, 0, sizeof(cl_mem), gaugefield->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(overrelax_odd, 2, sizeof(cl_mem), prng->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	for(cl_int i = 0; i < NDIM; i++) {
		clerr = clSetKernelArg(overrelax_odd, 1, sizeof(cl_int), &i);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
		get_device()->enqueue_kernel(overrelax_odd, global_work_size, ls);
	}
}

void hardware::code::Heatbath::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
	Opencl_Module::get_work_sizes(kernel, ls, gs, num_groups);

	//Query kernel name
	string kernelname = get_kernel_name(kernel);

	//Query specific sizes for kernels if needed
	//all of the following kernels are called with EnqueueKernel(gs), ls, num_groups are not needed!
	if (kernelname.compare("heatbath_even") == 0 || kernelname.compare("heatbath_odd") == 0 || kernelname.compare("overrelax_even") == 0 || kernelname.compare("overrelax_odd") == 0) {
		if( get_device()->get_device_type() == CL_DEVICE_TYPE_GPU ) {
			*gs = std::min(meta::get_volspace(get_parameters()) * get_parameters().get_ntime() / 2, hardware::buffers::get_prng_buffer_size(get_device()));
		} else {
			*gs = std::min(get_device()->get_num_compute_units(), hardware::buffers::get_prng_buffer_size(get_device()));
		}
		*ls = get_device()->get_preferred_local_thread_num();
		*num_groups = *gs / *ls;
	}
	return;
}

size_t hardware::code::Heatbath::get_read_write_size(const std::string& in) const
{
	//Depending on the compile-options, one has different sizes...
	size_t D = meta::get_float_size(get_parameters());
	size_t R = meta::get_mat_size(get_parameters());
	size_t S;
	//factor for complex numbers
	int C = 2;
	const size_t VOL4D = meta::get_vol4d(get_parameters());
	if(get_parameters().get_use_eo() == 1)
		S = meta::get_eoprec_spinorfieldsize(get_parameters());
	else
		S = meta::get_spinorfieldsize(get_parameters());
	//this is the same as in the function above
	if ( (in == "heatbath_even" ) || (in == "heatbath_odd") || (in == "overrelax_even") || (in == "overrelax_odd")) {
		//this kernel reads ingredients for 1 staple plus 1 su3matrix and writes 1 su3-matrix
		return VOL4D / 2 * C * D * R * (6 * (NDIM - 1) + 1 + 1 );
	}
	return 0;
}

uint64_t hardware::code::Heatbath::get_flop_size(const std::string& in) const
{
	const size_t VOL4D = meta::get_vol4d(get_parameters());
	uint64_t S;
	if(get_parameters().get_use_eo() == 1)
		S = meta::get_eoprec_spinorfieldsize(get_parameters());
	else
		S = meta::get_spinorfieldsize(get_parameters());
	//this is the same as in the function above
	///@NOTE: I do not distinguish between su3 and 3x3 matrices. This is a difference if one use e.g. REC12, but here one wants to have the "netto" flops for comparability.
	if ( (in == "heatbath_even" ) || (in == "heatbath_odd") ) {
		//this kernel calculates 1 staple (= 4*ND-1 su3_su3 + 2_ND-1 su3_add) plus NC*(2*su3_su3 80 flops for the su2 update)
		return VOL4D / 2 * (4 * (NDIM - 1) * meta::get_flop_su3_su3() + 2 * (NDIM - 1) * 18 + NC * (2 * meta::get_flop_su3_su3() + 80));
	}
	if ( (in == "overrelax_even") || (in == "overrelax_odd")) {
		//this kernel calculates 1 staple (= 4*ND-1 su3_su3 + 2_ND-1 su3_add) plus NC*(2*su3_su3 58 flops for the su2 update)
		return VOL4D / 2 * (4 * (NDIM - 1) * meta::get_flop_su3_su3() + 2 * (NDIM - 1) * 18 + NC * (2 * meta::get_flop_su3_su3() + 58));
	}
	return 0;
}

void hardware::code::Heatbath::print_profiling(const std::string& filename, int number) const
{
	Opencl_Module::print_profiling(filename, number);
	Opencl_Module::print_profiling(filename, heatbath_even);
	Opencl_Module::print_profiling(filename, heatbath_odd);
	Opencl_Module::print_profiling(filename, overrelax_even);
	Opencl_Module::print_profiling(filename, overrelax_odd);
}

hardware::code::Heatbath::Heatbath(const meta::Inputparameters& params, hardware::Device * device)
	: Opencl_Module(params, device)
{
	fill_kernels();
}

hardware::code::Heatbath::~Heatbath()
{
	clear_kernels();
}
