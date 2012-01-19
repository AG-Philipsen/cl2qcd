#include "opencl_module_spinors.h"

#include <algorithm>
#include <boost/regex.hpp>

#include "logger.hpp"

using namespace std;


void Opencl_Module_Spinors::fill_collect_options(stringstream* collect_options)
{
	Opencl_Module_Ran::fill_collect_options(collect_options);
	*collect_options << " -D_FERMIONS_" << " -DSPINORSIZE=" << get_parameters()->get_spinorsize() << " -DHALFSPINORSIZE=" << get_parameters()->get_halfspinorsize()
	                 << " -DSPINORFIELDSIZE=" << get_parameters()->get_spinorfieldsize() << " -DEOPREC_SPINORFIELDSIZE=" << get_parameters()->get_eoprec_spinorfieldsize();
	return;
}


void Opencl_Module_Spinors::fill_buffers()
{
	Opencl_Module_Ran::fill_buffers();

	//scratch buffers will be created on demand
	clmem_global_squarenorm_buf_glob = 0;
	clmem_scalar_product_buf_glob = 0;


	return;
}

void Opencl_Module_Spinors::fill_kernels()
{
	Opencl_Module_Ran::fill_kernels();

	basic_fermion_code = basic_opencl_code << "types_fermions.h" << "operations_su3vec.cl" << "operations_spinor.cl" << "spinorfield.cl";

	set_spinorfield_cold = createKernel("set_spinorfield_cold") << basic_fermion_code << "spinorfield_cold.cl";
	saxpy = createKernel("saxpy") << basic_fermion_code << "spinorfield_saxpy.cl";
	sax = createKernel("sax") << basic_fermion_code << "spinorfield_sax.cl";
	saxsbypz = createKernel("saxsbypz") << basic_fermion_code << "spinorfield_saxsbypz.cl";
	scalar_product = createKernel("scalar_product") << basic_fermion_code << "spinorfield_scalar_product.cl";
	scalar_product_reduction = createKernel("scalar_product_reduction") << basic_fermion_code << "spinorfield_scalar_product.cl";
	set_zero_spinorfield = createKernel("set_zero_spinorfield") << basic_fermion_code << "spinorfield_set_zero.cl";
	global_squarenorm = createKernel("global_squarenorm") << basic_fermion_code << "spinorfield_squarenorm.cl";
	global_squarenorm_reduction = createKernel("global_squarenorm_reduction") << basic_fermion_code << "spinorfield_squarenorm.cl";

	ratio = createKernel("ratio") << basic_opencl_code << "complex_ratio.cl";
	product = createKernel("product") << basic_opencl_code << "complex_product.cl";

	if(get_parameters()->get_use_eo() == true) {
		convert_from_eoprec = createKernel("convert_from_eoprec") << basic_fermion_code << "spinorfield_eo_convert.cl";
		convert_to_eoprec = createKernel("convert_to_eoprec") << basic_fermion_code << "spinorfield_eo_convert.cl";
		set_eoprec_spinorfield_cold = createKernel("set_eoprec_spinorfield_cold") << basic_fermion_code << "spinorfield_eo_cold.cl";
		saxpy_eoprec = createKernel("saxpy_eoprec") << basic_fermion_code << "spinorfield_eo_saxpy.cl";
		sax_eoprec = createKernel("sax_eoprec") << basic_fermion_code << "spinorfield_eo_sax.cl";
		saxsbypz_eoprec = createKernel("saxsbypz_eoprec") << basic_fermion_code << "spinorfield_eo_saxsbypz.cl";
		scalar_product_eoprec = createKernel("scalar_product_eoprec") << basic_fermion_code << "spinorfield_eo_scalar_product.cl";
		set_zero_spinorfield_eoprec = createKernel("set_zero_spinorfield_eoprec") << basic_fermion_code << "spinorfield_eo_zero.cl";
		global_squarenorm_eoprec = createKernel("global_squarenorm_eoprec") << basic_fermion_code << "spinorfield_eo_squarenorm.cl";
	}

	return;
}

void Opencl_Module_Spinors::clear_kernels()
{
	Opencl_Module_Ran::clear_kernels();

	cl_int clerr = CL_SUCCESS;

	clerr = clReleaseKernel(set_spinorfield_cold);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);

	clerr = clReleaseKernel(saxpy);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(sax);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(saxsbypz);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(scalar_product);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(scalar_product_reduction);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);

	clerr = clReleaseKernel(set_zero_spinorfield);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(global_squarenorm);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(global_squarenorm_reduction);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(ratio);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(product);

	if(get_parameters()->get_use_eo()) {
		clerr = clReleaseKernel(saxpy_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(sax_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(scalar_product_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(set_zero_spinorfield_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(global_squarenorm_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}



	return;
}

void Opencl_Module_Spinors::clear_buffers()
{
	Opencl_Module_Ran::clear_buffers();

	cl_uint clerr = CL_SUCCESS;
	if(clmem_scalar_product_buf_glob != 0 ) {
		clerr = clReleaseMemObject(clmem_scalar_product_buf_glob);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
	}
	if(clmem_global_squarenorm_buf_glob != 0) {
		clerr = clReleaseMemObject(clmem_global_squarenorm_buf_glob);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clMemObject", __FILE__, __LINE__);
	}

	return;
}


void Opencl_Module_Spinors::get_work_sizes(const cl_kernel kernel, cl_device_type dev_type, size_t * ls, size_t * gs, cl_uint * num_groups)
{
	Opencl_Module_Ran::get_work_sizes(kernel, dev_type, ls, gs, num_groups);

	//Query specific sizes for kernels if needed
	string kernelname = get_kernel_name(kernel);
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


	return;
}

void Opencl_Module_Spinors::convert_from_eoprec_device(cl_mem in1, cl_mem in2, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(convert_from_eoprec, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(convert_from_eoprec, 0, sizeof(cl_mem), &in1);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(convert_from_eoprec, 1, sizeof(cl_mem), &in2);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(convert_from_eoprec, 2, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(convert_from_eoprec , gs2, ls2);
}

void Opencl_Module_Spinors::convert_to_eoprec_device(cl_mem out1, cl_mem out2, cl_mem in)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(convert_to_eoprec, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(convert_to_eoprec, 0, sizeof(cl_mem), &out1);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(convert_to_eoprec, 1, sizeof(cl_mem), &out2);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(convert_to_eoprec, 2, sizeof(cl_mem), &in);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(convert_to_eoprec , gs2, ls2);
}

//BLAS-functions
void Opencl_Module_Spinors::saxpy_device(cl_mem x, cl_mem y, cl_mem alpha, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxpy, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxpy, 0, sizeof(cl_mem), &x);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy, 1, sizeof(cl_mem), &y);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy, 2, sizeof(cl_mem), &alpha);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy, 3, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(saxpy , gs2, ls2);
}

void Opencl_Module_Spinors::sax_device(cl_mem x, cl_mem alpha, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(sax, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(sax, 0, sizeof(cl_mem), &x);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(sax, 1, sizeof(cl_mem), &alpha);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(sax, 2, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(sax , gs2, ls2);
}

void Opencl_Module_Spinors::set_spinorfield_cold_device(cl_mem inout)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(set_spinorfield_cold, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(set_spinorfield_cold, 0, sizeof(cl_mem), &inout);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(set_spinorfield_cold , gs2, ls2);
}

void Opencl_Module_Spinors::set_eoprec_spinorfield_cold_device(cl_mem inout)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(set_eoprec_spinorfield_cold, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(set_eoprec_spinorfield_cold, 0, sizeof(cl_mem), &inout);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(set_eoprec_spinorfield_cold , gs2, ls2);
}

void Opencl_Module_Spinors::saxpy_eoprec_device(cl_mem x, cl_mem y, cl_mem alpha, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxpy_eoprec, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxpy_eoprec, 0, sizeof(cl_mem), &x);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_eoprec, 1, sizeof(cl_mem), &y);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_eoprec, 2, sizeof(cl_mem), &alpha);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_eoprec, 3, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( saxpy_eoprec, gs2, ls2);
}

void Opencl_Module_Spinors::sax_eoprec_device(cl_mem x, cl_mem alpha, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(sax_eoprec, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(sax_eoprec, 0, sizeof(cl_mem), &x);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(sax_eoprec, 1, sizeof(cl_mem), &alpha);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(sax_eoprec, 2, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( sax_eoprec, gs2, ls2);
}

void Opencl_Module_Spinors::saxsbypz_device(cl_mem x, cl_mem y, cl_mem z, cl_mem alpha, cl_mem beta, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxsbypz, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxsbypz, 0, sizeof(cl_mem), &x);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxsbypz, 1, sizeof(cl_mem), &y);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxsbypz, 2, sizeof(cl_mem), &z);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxsbypz, 3, sizeof(cl_mem), &alpha);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxsbypz, 4, sizeof(cl_mem), &beta);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxsbypz, 5, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( saxsbypz, gs2, ls2);
}

void Opencl_Module_Spinors::saxsbypz_eoprec_device(cl_mem x, cl_mem y, cl_mem z, cl_mem alpha, cl_mem beta, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxsbypz_eoprec, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxsbypz_eoprec, 0, sizeof(cl_mem), &x);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxsbypz_eoprec, 1, sizeof(cl_mem), &y);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxsbypz_eoprec, 2, sizeof(cl_mem), &z);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxsbypz_eoprec, 3, sizeof(cl_mem), &alpha);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxsbypz_eoprec, 4, sizeof(cl_mem), &beta);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxsbypz_eoprec, 5, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( saxsbypz_eoprec, gs2, ls2);
}

void Opencl_Module_Spinors::set_complex_to_scalar_product_device(cl_mem a, cl_mem b, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(scalar_product, this->get_device_type(), &ls2, &gs2, &num_groups);
	int global_buf_size_complex = sizeof(hmc_complex) * num_groups;
	if( clmem_scalar_product_buf_glob == 0 ) clmem_scalar_product_buf_glob = create_rw_buffer(global_buf_size_complex);

	//set arguments
	int clerr = clSetKernelArg(scalar_product, 0, sizeof(cl_mem), &a);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product, 1, sizeof(cl_mem), &b);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product, 2, sizeof(cl_mem), &clmem_scalar_product_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product, 3, sizeof(hmc_complex) * ls2, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(scalar_product , gs2, ls2);

	/// @todo Here the wait is needed. Replace this call by a clWaitForEvents!
	clerr = clFinish(get_queue());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clFinish", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product_reduction, 0, sizeof(cl_mem), &clmem_scalar_product_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product_reduction, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( scalar_product_reduction, gs2, ls2);
}

void Opencl_Module_Spinors::set_complex_to_scalar_product_eoprec_device(cl_mem a, cl_mem b, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(scalar_product_eoprec, this->get_device_type(), &ls2, &gs2, &num_groups);
	int global_buf_size_complex = sizeof(hmc_complex) * num_groups;
	if( clmem_scalar_product_buf_glob == 0 ) clmem_scalar_product_buf_glob = create_rw_buffer(global_buf_size_complex);

	//set arguments
	int clerr = clSetKernelArg(scalar_product_eoprec, 0, sizeof(cl_mem), &a);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product_eoprec, 1, sizeof(cl_mem), &b);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product_eoprec, 2, sizeof(cl_mem), &clmem_scalar_product_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product_eoprec, 3, sizeof(hmc_complex) * ls2, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( scalar_product_eoprec, gs2, ls2);
	clerr = clFinish(get_queue());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clFinish", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product_reduction, 0, sizeof(cl_mem), &clmem_scalar_product_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product_reduction, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( scalar_product_reduction, gs2, ls2);
}


void Opencl_Module_Spinors::set_complex_to_ratio_device(cl_mem a, cl_mem b, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(ratio, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(ratio, 0, sizeof(cl_mem), &a);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(ratio, 1, sizeof(cl_mem), &b);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(ratio, 2, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( ratio, gs2, ls2);
}

void Opencl_Module_Spinors::set_complex_to_product_device(cl_mem a, cl_mem b, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(product, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(product, 0, sizeof(cl_mem), &a);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(product, 1, sizeof(cl_mem), &b);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(product, 2, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(product , gs2, ls2);
}

void Opencl_Module_Spinors::set_float_to_global_squarenorm_device(cl_mem a, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(global_squarenorm, this->get_device_type(), &ls2, &gs2, &num_groups);
	int global_buf_size_float = sizeof(hmc_float) * num_groups;
	if( clmem_global_squarenorm_buf_glob == 0 ) clmem_global_squarenorm_buf_glob = create_rw_buffer(global_buf_size_float);

	//set arguments
	int clerr = clSetKernelArg(global_squarenorm, 0, sizeof(cl_mem), &a);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	//CP: these do not have to be args of the function since they are global objects to the class opencl??
	clerr = clSetKernelArg(global_squarenorm, 1, sizeof(cl_mem), &clmem_global_squarenorm_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(global_squarenorm, 2, sizeof(hmc_float) * ls2, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(global_squarenorm , gs2, ls2);

	/// @todo Here the wait is needed. Replace this call by a clWaitForEvents!
	clerr = clFinish(get_queue());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clFinish", __FILE__, __LINE__);

	clerr = clSetKernelArg(global_squarenorm_reduction, 0, sizeof(cl_mem), &clmem_global_squarenorm_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(global_squarenorm_reduction, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( global_squarenorm_reduction, gs2, ls2);
}

void Opencl_Module_Spinors::set_float_to_global_squarenorm_eoprec_device(cl_mem a, cl_mem out)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(global_squarenorm_eoprec, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int global_buf_size_float = sizeof(hmc_float) * num_groups;
	if( clmem_global_squarenorm_buf_glob == 0 ) clmem_global_squarenorm_buf_glob = create_rw_buffer(global_buf_size_float);

	int clerr = clSetKernelArg(global_squarenorm_eoprec, 0, sizeof(cl_mem), &a);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	//CP: these do not have to be args of the function since they are global objects to the class opencl??
	clerr = clSetKernelArg(global_squarenorm_eoprec, 1, sizeof(cl_mem), &clmem_global_squarenorm_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(global_squarenorm_eoprec, 2, sizeof(hmc_float) * ls2, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( global_squarenorm_eoprec, gs2, ls2);

	/// @todo Here the wait is needed. Replace this call by a clWaitForEvents!
	clFinish(get_queue());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clFinish", __FILE__, __LINE__);

	clerr = clSetKernelArg(global_squarenorm_reduction, 0, sizeof(cl_mem), &clmem_global_squarenorm_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(global_squarenorm_reduction, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( global_squarenorm_reduction, gs2, ls2);
}

void Opencl_Module_Spinors::set_zero_spinorfield_device(cl_mem x)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(set_zero_spinorfield, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(set_zero_spinorfield, 0, sizeof(cl_mem), &x);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(set_zero_spinorfield , gs2, ls2);
}

void Opencl_Module_Spinors::set_zero_spinorfield_eoprec_device(cl_mem x)
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(set_zero_spinorfield_eoprec, this->get_device_type(), &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(set_zero_spinorfield_eoprec, 0, sizeof(cl_mem), &x);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel( set_zero_spinorfield_eoprec, gs2, ls2);
}


#ifdef _PROFILING_
usetimer* Opencl_Module_Spinors::get_timer(const char * in)
{
	usetimer *noop = NULL;
	noop = Opencl_Module_Ran::get_timer(in);
	if(noop != NULL) return noop;
	if (strcmp(in, "set_spinorfield_cold") == 0) {
		return &this->timer_set_spinorfield_cold;
	}
	if (strcmp(in, "set_eoprec_spinorfield_cold") == 0) {
		return &this->timer_set_eoprec_spinorfield_cold;
	}
	if (strcmp(in, "convert_from_eoprec") == 0) {
		return &this->timer_convert_from_eoprec;
	}
	if (strcmp(in, "convert_to_eoprec") == 0) {
		return &this->timer_convert_to_eoprec;
	}
	if (strcmp(in, "saxpy") == 0) {
		return &(this->timer_saxpy);
	}
	if (strcmp(in, "sax") == 0) {
		return &(this->timer_sax);
	}
	if (strcmp(in, "saxsbypz") == 0) {
		return &this->timer_saxsbypz;
	}
	if (strcmp(in, "set_zero_spinorfield") == 0) {
		return &this->timer_set_zero_spinorfield;
	}
	if (strcmp(in, "saxpy_eoprec") == 0) {
		return &this->timer_saxpy_eoprec;
	}
	if (strcmp(in, "sax_eoprec") == 0) {
		return &this->timer_sax_eoprec;
	}
	if (strcmp(in, "saxsbypz_eoprec") == 0) {
		return &this->timer_saxsbypz_eoprec;
	}
	if (strcmp(in, "set_zero_spinorfield_eoprec") == 0) {
		return &this->timer_set_zero_spinorfield_eoprec;
	}
	if (strcmp(in, "scalar_product") == 0) {
		return &this->timer_scalar_product;
	}
	if (strcmp(in, "scalar_product_reduction") == 0) {
		return &this->timer_scalar_product_reduction;
	}
	if (strcmp(in, "global_squarenorm") == 0) {
		return &this->timer_global_squarenorm;
	}
	if (strcmp(in, "global_squarenorm_reduction") == 0) {
		return &this->timer_global_squarenorm_reduction;
	}
	if (strcmp(in, "scalar_product_eoprec") == 0) {
		return &this->timer_scalar_product_eoprec;
	}
	if (strcmp(in, "global_squarenorm_eoprec") == 0) {
		return &this->timer_global_squarenorm_eoprec;
	}
	if (strcmp(in, "ratio") == 0) {
		return &this->timer_ratio;
	}
	if (strcmp(in, "product") == 0) {
		return &this->timer_product;
	} else {
		return NULL;
	}
}

int Opencl_Module_Spinors::get_read_write_size(const char * in, inputparameters * parameters)
{
	int result = Opencl_Module_Ran::get_read_write_size(in, parameters);
	if (result != 0) return result;
	//Depending on the compile-options, one has different sizes...
	int D = (*parameters).get_float_size();
	//this returns the number of entries in an su3-matrix
	int R = (*parameters).get_mat_size();
	int S = get_parameters()->get_spinorfieldsize();
	int Seo = get_parameters()->get_eoprec_spinorfieldsize();
	//factor for complex numbers
	int C = 2;
	//this is the same as in the function above
	//NOTE: 1 spinor has NC*NDIM = 12 complex entries
	if (strcmp(in, "set_spinorfield_cold") == 0) {
		//this kernel writes 1 spinor
		return C * 12 * D * S;
	}
	if (strcmp(in, "set_eoprec_spinorfield_cold") == 0) {
		//this kernel writes 1 spinor
		return C * 12 * D * Seo;
	}
	if (strcmp(in, "convert_from_eoprec") == 0) {
		//this kernel reads 2 spinor and writes 2 spinors per site
		///@todo is this right??
		return 2 * 2 * C * 12 * D * Seo;
	}
	if (strcmp(in, "convert_to_eoprec") == 0) {
		//this kernel reads 2 spinor and writes 2 spinors per site
		return 2 * 2 * C * 12 * D * Seo;
	}
	if (strcmp(in, "saxpy") == 0) {
		//this kernel reads 2 spinor, 2 complex number and writes 1 spinor per site
		return C * D * S * (12 * (2 + 1) + 2);
	}
	if (strcmp(in, "sax") == 0) {
		//this kernel reads 1 spinor, 1 complex number and writes 1 spinor per site
		return C * D * S * (12 * (1 + 1) + 1);
	}
	if (strcmp(in, "saxsbypz") == 0) {
		//this kernel reads 3 spinor, 2 complex number and writes 1 spinor per site
		return C * D * S * (12 * (3 + 1) + 2);
	}
	if (strcmp(in, "set_zero_spinorfield") == 0) {
		//this kernel writes 1 spinor
		return C * 12 * D * S;
	}
	if (strcmp(in, "saxpy_eoprec") == 0) {
		//this kernel reads 2 spinor, 1 complex number and writes 1 spinor per site
		return C * D * Seo * (12 * (2 + 1) + 1);
	}
	if (strcmp(in, "sax_eoprec") == 0) {
		//this kernel reads 1 spinor, 1 complex number and writes 1 spinor per site
		return C * D * Seo * (12 * (1 + 1) + 1);
	}
	if (strcmp(in, "saxsbypz_eoprec") == 0) {
		//this kernel reads 3 spinor, 2 complex number and writes 1 spinor per site
		return C * D * Seo * (12 * (3 + 1) + 2);
	}
	if (strcmp(in, "set_zero_spinorfield_eoprec") == 0) {
		//this kernel writes 1 spinor
		return C * 12 * D * Seo;
	}
	if (strcmp(in, "scalar_product") == 0) {
		//this kernel reads 2 spinors and writes 1 complex number
		/// @NOTE: here, the local reduction is not taken into account
		return C * D * S * ( 2 * 12  + 1 );
	}
	if (strcmp(in, "scalar_product_reduction") == 0) {
		//this kernel reads NUM_GROUPS complex numbers and writes 1 complex number
		//query work-sizes for kernel to get num_groups
		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(scalar_product_reduction, this->get_device_type(), &ls2, &gs2, &num_groups);
		return C * D * (num_groups + 1);
	}
	if (strcmp(in, "global_squarenorm") == 0) {
		//this kernel reads 1 spinor and writes 1 real number
		/// @NOTE: here, the local reduction is not taken into account
		return D * S * (C * 12  + 1 );
	}
	if (strcmp(in, "global_squarenorm_reduction") == 0) {
		//this kernel reads NUM_GROUPS real numbers and writes 1 real number
		//query work-sizes for kernel to get num_groups
		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(scalar_product_reduction, this->get_device_type(), &ls2, &gs2, &num_groups);
		return D * (num_groups + 1);
	}
	if (strcmp(in, "scalar_product_eoprec") == 0) {
		//this kernel reads 2 spinors and writes 1 complex number
		/// @NOTE: here, the local reduction is not taken into account
		return C * D * Seo * ( 2 * 12  + 1 );
	}
	if (strcmp(in, "global_squarenorm_eoprec") == 0) {
		//this kernel reads 1 spinor and writes 1 real number
		/// @NOTE: here, the local reduction is not taken into account
		return D * Seo * (C * 12  + 1 );
	}
	if (strcmp(in, "ratio") == 0) {
		//this kernel reads 2 complex numbers and writes 1 complex number
		return C * D * (2 + 1);
	}
	if (strcmp(in, "product") == 0) {
		//this kernel reads 2 complex numbers and writes 1 complex number
		return C * D * (2 + 1);
	}
	return 0;
}

int Opencl_Module_Spinors::get_flop_size(const char * in, inputparameters * parameters)
{
	int result = Opencl_Module_Ran::get_flop_size(in, parameters);
	if (result != 0) return result;
	int S = get_parameters()->get_spinorfieldsize();
	int Seo = get_parameters()->get_eoprec_spinorfieldsize();
	//this is the same as in the function above
	if (strcmp(in, "set_spinorfield_cold") == 0) {
		//this kernel performs 1. / sqrt((12.f * VOL4D)) and real_multiply_spinor for each site
		return S * ( 3 + 24);
	}
	if (strcmp(in, "set_eoprec_spinorfield_cold") == 0) {
		//this kernel performs 1. / sqrt((12.f * VOL4D/2)) and real_multiply_spinor for each site
		return Seo * ( 3 + 24);
	}
	if (strcmp(in, "convert_from_eoprec") == 0) {
		//this kernel does not perform any flop, he just copies memory
		return 0;
	}
	if (strcmp(in, "convert_to_eoprec") == 0) {
		//this kernel does not perform any flop, he just copies memory
		return 0;
	}
	if (strcmp(in, "saxpy") == 0) {
		//this kernel performs on each site spinor_times_complex and spinor_add
		return S * (NDIM * NC * ( get_parameters()->get_flop_complex_mult() + 2) );
	}
	if (strcmp(in, "sax") == 0) {
		//this kernel performs on each site spinor_times_complex
		return S * (NDIM * NC * ( get_parameters()->get_flop_complex_mult() ) );
	}
	if (strcmp(in, "saxsbypz") == 0) {
		//this kernel performs on each 2 * site spinor_times_complex and 2 * spinor_add
		return S * (NDIM * NC * 2 * ( get_parameters()->get_flop_complex_mult() + 2) );
	}
	if (strcmp(in, "set_zero_spinorfield") == 0) {
		//this kernel does not do any flop
		return 0;
	}
	if (strcmp(in, "saxpy_eoprec") == 0) {
		//this kernel performs on each site spinor_times_complex and spinor_add
		return Seo * (NDIM * NC * ( get_parameters()->get_flop_complex_mult() + 2) );
	}
	if (strcmp(in, "sax_eoprec") == 0) {
		//this kernel performs on each site spinor_times_complex
		return S * (NDIM * NC * ( get_parameters()->get_flop_complex_mult() ) );
	}
	if (strcmp(in, "saxsbypz_eoprec") == 0) {
		//this kernel performs on each 2 * site spinor_times_complex and 2 * spinor_add
		return Seo * (NDIM * NC * 2 * ( get_parameters()->get_flop_complex_mult() + 2) );
	}
	if (strcmp(in, "set_zero_spinorfield_eoprec") == 0) {
		//this kernel does not do any flop
		return 0;
	}
	if (strcmp(in, "scalar_product") == 0) {
		//this kernel performs spinor*spinor on each site and then adds S-1 complex numbers
		return S * get_parameters()->get_flop_spinor_spinor() + (S - 1) * 2;
	}
	if (strcmp(in, "scalar_product_reduction") == 0) {
		return 1000000000000000000000000;
	}
	if (strcmp(in, "global_squarenorm") == 0) {
		//this kernel performs spinor_squarenorm on each site and then adds S-1 complex numbers
		return S * get_parameters()->get_flop_spinor_sqnorm() + (S - 1) * 2;
	}
	if (strcmp(in, "global_squarenorm_reduction") == 0) {
		return 1000000000000000000000000;
	}
	if (strcmp(in, "scalar_product_eoprec") == 0) {
		//this kernel performs spinor*spinor on each site and then adds S-1 complex numbers
		return Seo * get_parameters()->get_flop_spinor_spinor() + (Seo - 1) * 2;
	}
	if (strcmp(in, "global_squarenorm_eoprec") == 0) {
		//this kernel performs spinor_squarenorm on each site and then adds S-1 complex numbers
		return Seo * get_parameters()->get_flop_spinor_sqnorm() + (Seo - 1) * 2;
	}
	if (strcmp(in, "ratio") == 0) {
		return 11;
	}
	if (strcmp(in, "product") == 0) {
		return get_parameters()->get_flop_complex_mult();
	}
	return 0;
}

void Opencl_Module_Spinors::print_profiling(std::string filename, int number)
{
	Opencl_Module_Ran::print_profiling(filename, number);
	const char * kernelName;
	kernelName = "set_spinorfield_cold";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "set_eoprec_spinorfield_cold";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "convert_from_eoprec";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "convert_to_eoprec";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "saxpy";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "sax";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "saxsbypz";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "set_zero_spinorfield";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "saxpy_eoprec";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "sax_eoprec";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "saxsbypz_eoprec";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "set_zero_spinorfield_eoprec";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "scalar_product";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "scalar_product_reduction";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "global_squarenorm";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "global_squarenorm_reduction";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "scalar_product_eoprec";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "global_squarenorm_eoprec";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "ratio";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );
	kernelName = "product";
	Opencl_Module_Ran::print_profiling(filename, kernelName, (*this->get_timer(kernelName)).getTime(), (*this->get_timer(kernelName)).getNumMeas(), this->get_read_write_size(kernelName, parameters), this->get_flop_size(kernelName, parameters) );

}
#endif
