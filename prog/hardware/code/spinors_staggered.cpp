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
	if(check_su3vec_for_SOA(device)) {
		options << " -D EOPREC_SPINORFIELD_STRIDE=" << get_su3vec_buffer_stride(get_eoprec_spinorfieldsize(mem_size), device);
	}

	return options.str();
}

void hardware::code::Spinors_staggered::fill_kernels()
{
	basic_fermion_code = get_device()->get_gaugefield_code()->get_sources() << ClSourcePackage(collect_build_options(get_device(), get_parameters())) << "types_fermions.h" << "operations_su3vec.cl" << "operations_spinor.cl" << "spinorfield_staggered.cl";
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
	set_gaussian_spinorfield_stagg = createKernel("set_gaussian_spinorfield_stagg") << basic_fermion_code << prng_code << "spinorfield_staggered_gaussian.cl";
	//Complex number operations
	convert_stagg = createKernel("convert_float_to_complex") << get_device()->get_gaugefield_code()->get_sources() << "complex_convert.cl";
	ratio_stagg = createKernel("ratio") << get_device()->get_gaugefield_code()->get_sources() << "complex_ratio.cl";
	product_stagg = createKernel("product") << get_device()->get_gaugefield_code()->get_sources() << "complex_product.cl";
	//Fields algebra operations
	sax_stagg = createKernel("sax_staggered") << basic_fermion_code << "spinorfield_staggered_sax.cl";
	saxpy_stagg = createKernel("saxpy_staggered") << basic_fermion_code << "spinorfield_staggered_saxpy.cl";
	saxpbypz_stagg = createKernel("saxpbypz_staggered") << basic_fermion_code << "spinorfield_staggered_saxpbypz.cl";
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
	clerr = clReleaseKernel(set_gaussian_spinorfield_stagg);
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
	clerr = clReleaseKernel(saxpy_stagg);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(saxpbypz_stagg);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
}


void hardware::code::Spinors_staggered::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
	Opencl_Module::get_work_sizes(kernel, ls, gs, num_groups);

	// kernels that use random numbers must not exceed the size of the random state array
	if(kernel == set_gaussian_spinorfield_stagg
	   /*|| kernel == generate_gaussian_spinorfield_eo*/) {
		if(*gs > hardware::buffers::get_prng_buffer_size(get_device(), get_parameters())) {
			*gs = hardware::buffers::get_prng_buffer_size(get_device(), get_parameters());
			logger.error() << "I changed gs without changing neither ls nor num_groups (in Spinors_staggered::get_work_sizes)!!!";
		}
	}
	
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


void hardware::code::Spinors_staggered::saxpy_device(const hardware::buffers::Plain<su3vec> * x, const hardware::buffers::Plain<su3vec> * y, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Plain<su3vec> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxpy_stagg, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxpy_stagg, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_stagg, 1, sizeof(cl_mem), y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_stagg, 2, sizeof(cl_mem), alpha->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_stagg, 3, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(saxpy_stagg , gs2, ls2);
}


void hardware::code::Spinors_staggered::saxpbypz_device(const hardware::buffers::Plain<su3vec> * x, const hardware::buffers::Plain<su3vec> * y, const hardware::buffers::Plain<su3vec> * z, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Plain<hmc_complex> * beta, const hardware::buffers::Plain<su3vec> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxpbypz_stagg, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxpbypz_stagg, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpbypz_stagg, 1, sizeof(cl_mem), y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpbypz_stagg, 2, sizeof(cl_mem), z->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpbypz_stagg, 3, sizeof(cl_mem), alpha->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpbypz_stagg, 4, sizeof(cl_mem), beta->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpbypz_stagg, 5, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(saxpbypz_stagg, gs2, ls2);
}


void hardware::code::Spinors_staggered::set_gaussian_spinorfield_device(const hardware::buffers::Plain<su3vec> * in, const hardware::buffers::PRNGBuffer * prng) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(set_gaussian_spinorfield_stagg, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(set_gaussian_spinorfield_stagg, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(set_gaussian_spinorfield_stagg, 1, sizeof(cl_mem), prng->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(set_gaussian_spinorfield_stagg  , gs2, ls2);

	if(logger.beDebug()) {
		hardware::buffers::Plain<hmc_float> force_tmp(1, get_device());
		hmc_float resid;
		get_device()->get_spinor_staggered_code()->set_float_to_global_squarenorm_device(in, &force_tmp);
		force_tmp.dump(&resid);
		logger.debug() <<  "\tinit gaussian spinorfield:\t" << resid;
		if(resid != resid) {
			throw Print_Error_Message("calculation of gaussian spinorfield gave nan! Aborting...", __FILE__, __LINE__);
		}
	}
}



/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////


size_t hardware::code::Spinors_staggered::get_read_write_size(const std::string& in) const
{
	//Depending on the compile-options, one has different sizes...
	size_t D = meta::get_float_size(get_parameters());
	size_t S = get_spinorfieldsize(get_parameters());
	size_t Seo = get_eoprec_spinorfieldsize(get_parameters());
	//factor for complex numbers
	int C = 2;
	//NOTE: 1 spinor has NC*NSPIN = 3*1 = 3 complex entries
	if (in == "global_squarenorm_staggered") {
		//this kernel reads 1 su3vec and writes 1 real number
		/// @NOTE: here, the local reduction is not taken into account
		return D * S * (C * NC + 1);
	}
	if (in == "global_squarenorm_reduction") {
		//this kernel reads NUM_GROUPS real numbers and writes 1 real number
		//query work-sizes for kernel to get num_groups
		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(global_squarenorm_stagg, &ls2, &gs2, &num_groups);
		return D * (num_groups + 1);
	}
	if (in == "scalar_product_staggered") {
		//this kernel reads 2 spinors and writes 1 complex number
		/// @NOTE: here, the local reduction is not taken into account
		return C * D * S * (NC * 2 + 1);
	}
	if (in == "scalar_product_reduction") {
		//this kernel reads NUM_GROUPS complex numbers and writes 1 complex number
		//query work-sizes for kernel to get num_groups
		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(scalar_product_stagg, &ls2, &gs2, &num_groups);
		return C * D * (num_groups + 1);
	}
	if (in == "set_zero_spinorfield_stagg") {
		//this kernel writes 1 su3vec
		return C * D * S * NC;
	}
	if (in == "set_cold_spinorfield_stagg") {
		//this kernel writes 1 su3vec
		return C * D * S * NC;
	}
	if (in == "convert_float_to_complex") {
		//this kernel reads 1 float and writes 1 complex number
		return (C + 1) * D;
	}
	if (in == "set_gaussian_spinorfield_stagg") {
		//this kernel writes 1 su3vec per site
		return ( NC * C ) * D * S;
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
		//this kernel reads 1 su3vec, 1 complex number and writes 1 su3vec per site
		return C * D * S * (NC * (1 + 1) + 1);
	}
	if (in == "saxpy_staggered") {
		//this kernel reads 2 su3vec, 2 complex number and writes 1 su3vec per site
		return C * D * S * (NC * (2 + 1) + 2);
	}
	if (in == "saxpbypz_staggered") {
		//this kernel reads 3 su3vec, 2 complex number and writes 1 su3vec per site
		return C * D * S * (NC * (3 + 1) + 2);
	}
	return 0;
}

uint64_t hardware::code::Spinors_staggered::get_flop_size(const std::string& in) const
{
	uint64_t S = get_spinorfieldsize(get_parameters());
	uint64_t Seo = get_eoprec_spinorfieldsize(get_parameters());
	//this is the same as in the function above
	if (in == "global_squarenorm_staggered") {
		//this kernel performs spinor_squarenorm on each site and then adds S-1 complex numbers
		//Note that the sum of S numbers can be done in several way: the least efficient way is
		//to add the first 2 numbers, then the third, then the fourth and so on, performing S-1
		//additions. This is not what is done in the code but it is a good estimation because
		//we know for sure that the code will be a bit faster (it is somehow a boundary estimate)
		return S * meta::get_flop_su3vec_sqnorm() + (S - 1) * 2;
	}
	if (in == "global_squarenorm_reduction") {
		//This if should not be entered since the sum of the site squarenorms
		//has already taken into account with the (S-1)*2 term in the previous if
		return 1000000000000000000000000;
	}
	if (in == "scalar_product_staggered") {
		//this kernel performs su3vec_scalarproduct on each site and then adds S-1 complex numbers
		//Note that the sum of S numbers can be done in several way: the least efficient way is
		//to add the first 2 numbers, then the third, then the fourth and so on, performing S-1
		//additions. This is not what is done in the code but it is a good estimation because
		//we know for sure that the code will be a bit faster (it is somehow a boundary estimate)
		return S * meta::get_flop_su3vec_su3vec() + (S - 1) * 2;
	}
	if (in == "scalar_product_reduction") {
		//This if should not be entered since the sum of the site squarenorms
		//has already taken into account with the (S-1)*2 term in the previous if
		return 1000000000000000000000000;
	}
	if (in == "set_zero_spinorfield_stagg") {
		//this kernel does not do any flop
		return 0;
	}
	if (in == "set_cold_spinorfield_stagg") {
		//this kernel performs 1. / sqrt((3.f * VOL4D)) and su3vec_times_real for each site
		return S * (3 + NC * 2);
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
	if (in == "set_gaussian_spinorfield_stagg") {
		//this kernel performs NC multiplications per site
		///@todo ? I did not count the gaussian normal pair production, which is very complicated...
		return NC * S;
	}
	if (in == "sax_staggered") {
		//this kernel performs on each site su3vec_times_complex
		return S * (NC * (meta::get_flop_complex_mult()));
	}
	if (in == "saxpy_staggered") {
		//this kernel performs on each site su3vec_times_complex and su3vec_acc
		return S * (NC * (meta::get_flop_complex_mult() + 2));
	}
	if (in == "saxpbypz_staggered") {
		//this kernel performs on each 2 * site su3vec_times_complex and 2 * su3vec_acc
		return S * (NC * 2 * ( meta::get_flop_complex_mult() + 2));
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
	Opencl_Module::print_profiling(filename, set_gaussian_spinorfield_stagg);
	Opencl_Module::print_profiling(filename, ratio_stagg);
	Opencl_Module::print_profiling(filename, convert_stagg);
	Opencl_Module::print_profiling(filename, product_stagg);
	Opencl_Module::print_profiling(filename, sax_stagg);
	Opencl_Module::print_profiling(filename, saxpy_stagg);
	Opencl_Module::print_profiling(filename, saxpbypz_stagg);
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

