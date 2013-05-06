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
	//Squarenorm
	global_squarenorm_stagg = createKernel("global_squarenorm_staggered") << basic_fermion_code << "spinorfield_staggered_squarenorm.cl";
	global_squarenorm_reduction_stagg = createKernel("global_squarenorm_reduction") << basic_fermion_code << "spinorfield_staggered_squarenorm.cl";
	//Scalar Product
	scalar_product_stagg = createKernel("scalar_product_staggered") << basic_fermion_code << "spinorfield_staggered_scalar_product.cl";
	scalar_product_reduction_stagg = createKernel("scalar_product_reduction") << basic_fermion_code << "spinorfield_staggered_scalar_product.cl";
	//Setting fields
	set_zero_spinorfield_stagg = createKernel("set_zero_spinorfield_stagg") << basic_fermion_code << "spinorfield_staggered_set_zero.cl";
	set_cold_spinorfield_stagg = createKernel("set_cold_spinorfield_stagg") << basic_fermion_code << "spinorfield_staggered_set_cold.cl";
	//Complex number operations
	convert_stagg = createKernel("convert_float_to_complex") << get_device()->get_gaugefield_code()->get_sources() << "complex_convert.cl";
	ratio_stagg = createKernel("ratio") << get_device()->get_gaugefield_code()->get_sources() << "complex_ratio.cl";
	product_stagg = createKernel("product") << get_device()->get_gaugefield_code()->get_sources() << "complex_product.cl";
	//Fields algebra operations
	sax_stagg = createKernel("sax_staggered") << basic_fermion_code << "spinorfield_staggered_sax.cl";
}

void hardware::code::Spinors_staggered::clear_kernels()
{
	cl_int clerr = CL_SUCCESS;
	//Squarenorm
	clerr = clReleaseKernel(global_squarenorm_stagg);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(global_squarenorm_reduction_stagg);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	//Scalar Product
	clerr = clReleaseKernel(scalar_product_stagg);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(scalar_product_reduction_stagg);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	//Setting fields
	clerr = clReleaseKernel(set_zero_spinorfield_stagg);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(set_cold_spinorfield_stagg);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	//Complex numbers operations
	clerr = clReleaseKernel(convert_stagg);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(ratio_stagg);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(product_stagg);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	//Fields algebra operations
	clerr = clReleaseKernel(sax_stagg);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
}


void hardware::code::Spinors_staggered::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
	Opencl_Module::get_work_sizes(kernel, ls, gs, num_groups);

	//Query specific sizes for kernels if needed
	string kernelname = get_kernel_name(kernel);
	if(kernelname.compare("convert_float_to_complex") == 0) {
		*ls = 1;
		*gs = 1;
		*num_groups = 1;
	}
	if(kernelname.compare("ratio") == 0) {
		*ls = 1;
		*gs = 1;
		*num_groups = 1;
	}
	if(kernelname.compare("product") == 0) {
		*ls = 1;
		*gs = 1;
		*num_groups = 1;
	}
	//Whenever ls id manually modified, it is crucial to modify num_groups accordingly!
	if(kernel == global_squarenorm_stagg || kernel == scalar_product_stagg
           /*|| kernel == global_squarenorm_eoprec || kernel == scalar_product_eoprec */) {
	  if(*ls > 64) {
	    *ls = 64;
	    *num_groups = (*gs)/(*ls);
	  }
	  return;
	}
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void hardware::code::Spinors_staggered::global_squarenorm_reduction(const hardware::buffers::Plain<hmc_float> * out, const hardware::buffers::Plain<hmc_float> * tmp_buf) const
{
	cl_int clerr = clSetKernelArg(global_squarenorm_reduction_stagg, 0, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(global_squarenorm_reduction_stagg, 1, sizeof(cl_mem), tmp_buf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	cl_uint elems = tmp_buf->get_elements();
	clerr = clSetKernelArg(global_squarenorm_reduction_stagg, 2, sizeof(cl_uint), &elems);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(global_squarenorm_reduction_stagg, 1, 1);
}

void hardware::code::Spinors_staggered::set_float_to_global_squarenorm_device(const hardware::buffers::Plain<su3vec> * a, const hardware::buffers::Plain<hmc_float> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(global_squarenorm_stagg, &ls2, &gs2, &num_groups);

	hardware::buffers::Plain<hmc_float> tmp(num_groups, get_device());

	//set arguments
	int clerr = clSetKernelArg(global_squarenorm_stagg, 0, sizeof(cl_mem), a->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	//CP: these do not have to be args of the function since they are global objects to the class opencl??
	clerr = clSetKernelArg(global_squarenorm_stagg, 1, sizeof(cl_mem), tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(global_squarenorm_stagg, 2, sizeof(hmc_float) * ls2, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(global_squarenorm_stagg, gs2, ls2);

	global_squarenorm_reduction(out, &tmp);
}


void hardware::code::Spinors_staggered::set_complex_to_scalar_product_device(const hardware::buffers::Plain<su3vec> * a, const hardware::buffers::Plain<su3vec> * b, const hardware::buffers::Plain<hmc_complex> * out) const
{
  //query work-sizes for kernel
  size_t ls2, gs2;
  cl_uint num_groups;
  this->get_work_sizes(scalar_product_stagg, &ls2, &gs2, &num_groups);

  hardware::buffers::Plain<hmc_complex> tmp(num_groups, get_device());

  //set arguments
  int clerr = clSetKernelArg(scalar_product_stagg, 0, sizeof(cl_mem), a->get_cl_buffer());
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

  clerr = clSetKernelArg(scalar_product_stagg, 1, sizeof(cl_mem), b->get_cl_buffer());
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

  clerr = clSetKernelArg(scalar_product_stagg, 2, sizeof(cl_mem), tmp);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

  clerr = clSetKernelArg(scalar_product_stagg, 3, sizeof(hmc_complex) * ls2, static_cast<void*>(nullptr));
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

  get_device()->enqueue_kernel(scalar_product_stagg , gs2, ls2);


  clerr = clSetKernelArg(scalar_product_reduction_stagg, 0, sizeof(cl_mem), tmp);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

  clerr = clSetKernelArg(scalar_product_reduction_stagg, 1, sizeof(cl_mem), out->get_cl_buffer());
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

  clerr = clSetKernelArg(scalar_product_reduction_stagg, 2, sizeof(cl_uint), &num_groups);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

  get_device()->enqueue_kernel(scalar_product_reduction_stagg, gs2, ls2);
}


void hardware::code::Spinors_staggered::set_zero_spinorfield_device(const hardware::buffers::Plain<su3vec> * x) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(set_zero_spinorfield_stagg, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(set_zero_spinorfield_stagg, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(set_zero_spinorfield_stagg, gs2, ls2);
}


void hardware::code::Spinors_staggered::set_cold_spinorfield_device(const hardware::buffers::Plain<su3vec> * x) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(set_cold_spinorfield_stagg, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(set_cold_spinorfield_stagg, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(set_cold_spinorfield_stagg, gs2, ls2);
}


void hardware::code::Spinors_staggered::set_complex_to_float_device(const hardware::buffers::Plain<hmc_float> * in, const hardware::buffers::Plain<hmc_complex> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(convert_stagg, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(convert_stagg, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(convert_stagg, 1, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(convert_stagg, gs2, ls2);
}


void hardware::code::Spinors_staggered::set_complex_to_ratio_device(const hardware::buffers::Plain<hmc_complex> * a, const hardware::buffers::Plain<hmc_complex> * b, const hardware::buffers::Plain<hmc_complex> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(ratio_stagg, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(ratio_stagg, 0, sizeof(cl_mem), a->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(ratio_stagg, 1, sizeof(cl_mem), b->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(ratio_stagg, 2, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(ratio_stagg, gs2, ls2);
}


void hardware::code::Spinors_staggered::set_complex_to_product_device(const hardware::buffers::Plain<hmc_complex> * a, const hardware::buffers::Plain<hmc_complex> * b, const hardware::buffers::Plain<hmc_complex> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(product_stagg, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(product_stagg, 0, sizeof(cl_mem), a->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(product_stagg, 1, sizeof(cl_mem), b->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(product_stagg, 2, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(product_stagg, gs2, ls2);
}


void hardware::code::Spinors_staggered::sax_device(const hardware::buffers::Plain<su3vec> * x, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Plain<su3vec> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(sax_stagg, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(sax_stagg, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(sax_stagg, 1, sizeof(cl_mem), alpha->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(sax_stagg, 2, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(sax_stagg, gs2, ls2);
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
	if (in == "global_squarenorm_staggered") {
		return 1000000000000000000000000;
	}
	if (in == "global_squarenorm_reduction") {
		return 1000000000000000000000000;
	}
	if (in == "scalar_product_staggered") {
		return 1000000000000000000000000;
	}
	if (in == "scalar_product_reduction") {
		return 1000000000000000000000000;
	}
	if (in == "set_zero_spinorfield_stagg") {
		return 1000000000000000000000000;
	}
	if (in == "set_cold_spinorfield_stagg") {
		return 1000000000000000000000000;
	}
	if (in == "convert_float_to_complex") {
		//this kernel reads 1 float and writes 1 complex number
		return (C + 1) * D;
	}
	if (in == "ratio") {
		//this kernel reads 2 complex numbers and writes 1 complex number
		return C * D * (2 + 1);
	}
	if (in == "product") {
		//this kernel reads 2 complex numbers and writes 1 complex number
		return C * D * (2 + 1);
	}
	if (in == "sax_staggered") {
		return 1000000000000000000000000;
	}
	return 0;
}

uint64_t hardware::code::Spinors_staggered::get_flop_size(const std::string& in) const
{
	uint64_t S = get_spinorfieldsize(get_parameters());
	uint64_t Seo = get_eoprec_spinorfieldsize(get_parameters());
	//this is the same as in the function above
	if (in == "global_squarenorm_staggered") {
		return 1000000000000000000000000;
	}
	if (in == "global_squarenorm_reduction") {
		return 1000000000000000000000000;
	}
	if (in == "scalar_product_staggered") {
		return 1000000000000000000000000;
	}
	if (in == "scalar_product_reduction") {
		return 1000000000000000000000000;
	}
	if (in == "set_zero_spinorfield_stagg") {
		return 1000000000000000000000000;
	}
	if (in == "set_cold_spinorfield_stagg") {
		return 1000000000000000000000000;
	}
	if (in == "convert_float_to_complex") {
		return 0;
	}
	if (in == "ratio") {
		return 11;
	}
	if (in == "product") {
		return meta::get_flop_complex_mult();
	}
	if (in == "sax_staggered") {
		return 1000000000000000000000000;
	}
	return 0;
}

void hardware::code::Spinors_staggered::print_profiling(const std::string& filename, int number) const
{
	Opencl_Module::print_profiling(filename, number);
	Opencl_Module::print_profiling(filename, global_squarenorm_stagg);
	Opencl_Module::print_profiling(filename, global_squarenorm_reduction_stagg);
	Opencl_Module::print_profiling(filename, scalar_product_stagg);
	Opencl_Module::print_profiling(filename, scalar_product_reduction_stagg);
	Opencl_Module::print_profiling(filename, set_zero_spinorfield_stagg);
	Opencl_Module::print_profiling(filename, set_cold_spinorfield_stagg);
	Opencl_Module::print_profiling(filename, ratio_stagg);
	Opencl_Module::print_profiling(filename, convert_stagg);
	Opencl_Module::print_profiling(filename, product_stagg);
	Opencl_Module::print_profiling(filename, sax_stagg);
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

