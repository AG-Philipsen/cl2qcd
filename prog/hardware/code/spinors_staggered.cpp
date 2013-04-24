#include "spinors_staggered.hpp"

#include "../../logger.hpp"
#include "../../meta/util.hpp"
#include "../device.hpp"
#include <cassert>
#include "gaugefield.hpp"
#include "prng.hpp"
#include "spinors.hpp"

using namespace std;

static std::string collect_build_options(hardware::Device * device, const meta::Inputparameters& params);

static std::string collect_build_options(hardware::Device * device, const meta::Inputparameters& params)
{
	using namespace hardware::buffers;
	using namespace hardware::code;

	const size_4 mem_size = device->get_mem_lattice_size();
	const size_4 local_size = device->get_local_lattice_size();

	std::ostringstream options;
	options.precision(16);
	options << "-D _FERMIONS_";
	options << " -D SPINORFIELDSIZE_GLOBAL=" << get_spinorfieldsize(params) << " -D EOPREC_SPINORFIELDSIZE_GLOBAL=" << get_eoprec_spinorfieldsize(params);
	options << " -D SPINORFIELDSIZE_LOCAL=" << get_spinorfieldsize(local_size) << " -D EOPREC_SPINORFIELDSIZE_MEM=" << get_eoprec_spinorfieldsize(local_size);
	options << " -D SPINORFIELDSIZE_MEM=" << get_spinorfieldsize(mem_size) << " -D EOPREC_SPINORFIELDSIZE_MEM=" << get_eoprec_spinorfieldsize(mem_size);
	if(check_Spinor_for_SOA(device)) {
		options << " -D EOPREC_SPINORFIELD_STRIDE=" << get_Spinor_buffer_stride(get_eoprec_spinorfieldsize(mem_size), device);
	}

	return options.str();
}

void hardware::code::Spinors_staggered::fill_kernels()
{
	basic_fermion_code = get_device()->get_gaugefield_code()->get_sources() << ClSourcePackage(collect_build_options(get_device(), get_parameters())) << "types_fermions.h" << "operations_su3vec.cl" << "operations_spinor.cl" << "spinorfield.cl";
	ClSourcePackage prng_code = get_device()->get_prng_code()->get_sources();
	global_squarenorm = createKernel("global_squarenorm_staggered") << basic_fermion_code << "spinorfield_staggered_squarenorm.cl";
	_global_squarenorm_reduction = createKernel("global_squarenorm_reduction") << basic_fermion_code << "spinorfield_staggered_squarenorm.cl";
}

void hardware::code::Spinors_staggered::clear_kernels()
{
	cl_int clerr = CL_SUCCESS;

	clerr = clReleaseKernel(global_squarenorm);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(_global_squarenorm_reduction);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
}


void hardware::code::Spinors_staggered::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
	Opencl_Module::get_work_sizes(kernel, ls, gs, num_groups);

	if(kernel == global_squarenorm) {
		if(*ls > 64) {
			*ls = 64;
		}
		return;
	}
}


void hardware::code::Spinors_staggered::global_squarenorm_reduction(const hardware::buffers::Plain<hmc_float> * out, const hardware::buffers::Plain<hmc_float> * tmp_buf) const
{
	cl_int clerr = clSetKernelArg(_global_squarenorm_reduction, 0, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(_global_squarenorm_reduction, 1, sizeof(cl_mem), tmp_buf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	cl_uint elems = tmp_buf->get_elements();
	clerr = clSetKernelArg(_global_squarenorm_reduction, 2, sizeof(cl_uint), &elems);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(_global_squarenorm_reduction, 1, 1);
}

void hardware::code::Spinors_staggered::set_float_to_global_squarenorm_device(const hardware::buffers::Plain<su3vec> * a, const hardware::buffers::Plain<hmc_float> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(global_squarenorm, &ls2, &gs2, &num_groups);

	hardware::buffers::Plain<hmc_float> tmp(num_groups, get_device());

	//set arguments
	int clerr = clSetKernelArg(global_squarenorm, 0, sizeof(cl_mem), a->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	//CP: these do not have to be args of the function since they are global objects to the class opencl??
	clerr = clSetKernelArg(global_squarenorm, 1, sizeof(cl_mem), tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(global_squarenorm, 2, sizeof(hmc_float) * ls2, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(global_squarenorm , gs2, ls2);

	global_squarenorm_reduction(out, &tmp);
}


size_t hardware::code::Spinors_staggered::get_read_write_size(const std::string& in) const
{
	//Depending on the compile-options, one has different sizes...
	size_t D = meta::get_float_size(get_parameters());
	size_t S = get_spinorfieldsize(get_parameters());
	size_t Seo = get_eoprec_spinorfieldsize(get_parameters());
	//factor for complex numbers
	int C = 2;
	//this is the same as in the function above
	if (in == "global_squarenorm") {
		return 1000000000000000000000000;
	}
	if (in == "global_squarenorm_reduction") {
		return 1000000000000000000000000;
	}
	return 0;
}

uint64_t hardware::code::Spinors_staggered::get_flop_size(const std::string& in) const
{
	uint64_t S = get_spinorfieldsize(get_parameters());
	uint64_t Seo = get_eoprec_spinorfieldsize(get_parameters());
	//this is the same as in the function above
	if (in == "global_squarenorm") {
		return 1000000000000000000000000;
	}
	if (in == "global_squarenorm_reduction") {
		return 1000000000000000000000000;
	}
	return 0;
}

void hardware::code::Spinors_staggered::print_profiling(const std::string& filename, int number) const
{
	Opencl_Module::print_profiling(filename, number);
	Opencl_Module::print_profiling(filename, global_squarenorm);
	Opencl_Module::print_profiling(filename, _global_squarenorm_reduction);
}

hardware::code::Spinors_staggered::Spinors_staggered(const meta::Inputparameters& params, hardware::Device * device)
	: Opencl_Module(params, device)
{
	fill_kernels();
}

hardware::code::Spinors_staggered::~Spinors_staggered()
{
	clear_kernels();
}

ClSourcePackage hardware::code::Spinors_staggered::get_sources() const noexcept
{
	return basic_fermion_code;
}

