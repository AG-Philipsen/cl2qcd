/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "spinors.hpp"

#include "../../host_functionality/logger.hpp"
#include "../device.hpp"
#include <cassert>

#include "flopUtilities.hpp"
#include "flopUtilities.hpp"
#include "gaugefield.hpp"
#include "prng.hpp"

using namespace std;

void hardware::code::Spinors::fill_kernels()
{
	basic_fermion_code = get_basic_sources() <<  "operations_geometry.cl" << "operations_complex.h" << "types_fermions.h" << "operations_su3vec.cl" << "operations_spinor.cl" << "spinorfield.cl";
	if(kernelParameters->getUseEo()) {
		basic_fermion_code = basic_fermion_code << "operations_spinorfield_eo.cl";
	}
	
	ClSourcePackage prng_code = get_device()->getPrngCode()->get_sources();
	
	logger.debug() << "Creating Spinors kernels...";
	
	//Reductions are really small kernels, so few needed options loaded by hands
	_global_squarenorm_reduction = createKernel("global_squarenorm_reduction")  << ClSourcePackage("-I " + std::string(SOURCEDIR) + " -D _INKERNEL_" + ((kernelParameters->getPrecision() == 64) ? (std::string(" -D _USEDOUBLEPREC_") + " -D _DEVICE_DOUBLE_EXTENSION_KHR_") : "")) << "types.h" << "spinorfield_squarenorm_reduction.cl";
	scalar_product_reduction = createKernel("scalar_product_reduction") << ClSourcePackage("-I " + std::string(SOURCEDIR) + " -D _INKERNEL_" + ((kernelParameters->getPrecision() == 64) ? (std::string(" -D _USEDOUBLEPREC_") + " -D _DEVICE_DOUBLE_EXTENSION_KHR_") : "")) << "types.h" << "operations_complex.h" << "spinorfield_scalar_product_reduction.cl";
	
	if(kernelParameters->getUseEo() ) {
		generate_gaussian_spinorfield_eo = createKernel("generate_gaussian_spinorfield_eo") << basic_fermion_code << prng_code << "spinorfield_eo_gaussian.cl";
		convert_from_eoprec = createKernel("convert_from_eoprec") << basic_fermion_code << "spinorfield_eo_convert.cl";
		convert_to_eoprec = createKernel("convert_to_eoprec") << basic_fermion_code << "spinorfield_eo_convert.cl";
		set_eoprec_spinorfield_cold = createKernel("set_eoprec_spinorfield_cold") << basic_fermion_code << "spinorfield_eo_cold.cl";
		saxpy_eoprec = createKernel("saxpy_eoprec") << basic_fermion_code << "spinorfield_eo_saxpy.cl";
		saxpy_arg_eoprec = createKernel("saxpy_arg_eoprec") << basic_fermion_code << "spinorfield_eo_saxpy.cl";
		sax_eoprec = createKernel("sax_eoprec") << basic_fermion_code << "spinorfield_eo_sax.cl";
		saxsbypz_eoprec = createKernel("saxsbypz_eoprec") << basic_fermion_code << "spinorfield_eo_saxsbypz.cl";
		scalar_product_eoprec = createKernel("scalar_product_eoprec") << basic_fermion_code << "spinorfield_eo_scalar_product.cl";
		set_zero_spinorfield_eoprec = createKernel("set_zero_spinorfield_eoprec") << basic_fermion_code << "spinorfield_eo_zero.cl";
		global_squarenorm_eoprec = createKernel("global_squarenorm_eoprec") << basic_fermion_code << "spinorfield_eo_squarenorm.cl";
		convertSpinorfieldToSOA_eo = createKernel("convertSpinorfieldToSOA_eo") << basic_fermion_code << "spinorfield_eo_convert.cl";
		convertSpinorfieldFromSOA_eo = createKernel("convertSpinorfieldFromSOA_eo") << basic_fermion_code << "spinorfield_eo_convert.cl";
		//merged kernels
		if (kernelParameters->getUseMergeKernelsSpinor() == true) {
			saxpy_AND_squarenorm_eo = createKernel("saxpy_AND_squarenorm_eo") << basic_fermion_code << "spinorfield_eo_saxpy_AND_squarenorm.cl";
		} else {
			saxpy_AND_squarenorm_eo = 0;
		}
	} else {
		generate_gaussian_spinorfield_eo = 0;
		convert_from_eoprec = 0;
		convert_to_eoprec =0;
		set_eoprec_spinorfield_cold = 0;
		sax_eoprec = 0;
		saxpy_eoprec = 0;
		saxpy_arg_eoprec = 0;
		saxsbypz_eoprec = 0;
		scalar_product_eoprec = 0;
		set_zero_spinorfield_eoprec = 0;
		global_squarenorm_eoprec = 0;
		convertSpinorfieldToSOA_eo = 0;
		convertSpinorfieldFromSOA_eo = 0;
	}
	//Always build non eo-prec kernels
	generate_gaussian_spinorfield = createKernel("generate_gaussian_spinorfield") << basic_fermion_code << prng_code << "spinorfield_gaussian.cl";
	set_spinorfield_cold = createKernel("set_spinorfield_cold") << basic_fermion_code << "spinorfield_cold.cl";
	saxpy = createKernel("saxpy") << basic_fermion_code << "spinorfield_saxpy.cl";
	saxpy_arg = createKernel("saxpy_arg") << basic_fermion_code << "spinorfield_saxpy.cl";
	sax = createKernel("sax") << basic_fermion_code << "spinorfield_sax.cl";
	saxsbypz = createKernel("saxsbypz") << basic_fermion_code << "spinorfield_saxsbypz.cl";
	scalar_product = createKernel("scalar_product") << basic_fermion_code << "spinorfield_scalar_product.cl";
	set_zero_spinorfield = createKernel("set_zero_spinorfield") << basic_fermion_code << "spinorfield_set_zero.cl";
	global_squarenorm = createKernel("global_squarenorm") << basic_fermion_code << "spinorfield_squarenorm.cl";
}

void hardware::code::Spinors::clear_kernels()
{
	cl_int clerr = CL_SUCCESS;
	
	logger.debug() << "Clearing Spinors kernels...";

	//Reductions
	clerr = clReleaseKernel(_global_squarenorm_reduction);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(scalar_product_reduction);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);

	if(kernelParameters->getUseEo()) {
		clerr = clReleaseKernel(generate_gaussian_spinorfield_eo);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(saxpy_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(saxpy_arg_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(sax_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(scalar_product_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(set_zero_spinorfield_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(global_squarenorm_eoprec);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(convertSpinorfieldToSOA_eo);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
		clerr = clReleaseKernel(convertSpinorfieldFromSOA_eo);
		if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	}
	//Always build non eo-prec kernels
	clerr = clReleaseKernel(generate_gaussian_spinorfield);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(set_spinorfield_cold);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(saxpy);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(saxpy_arg);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(sax);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(saxsbypz);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(scalar_product);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(set_zero_spinorfield);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(global_squarenorm);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
}


void hardware::code::Spinors::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
	Opencl_Module::get_work_sizes(kernel, ls, gs, num_groups);

	// kernels that use random numbers must not exceed the size of the random state array
	if(kernel == generate_gaussian_spinorfield
	   || kernel == generate_gaussian_spinorfield_eo) {
		if(*gs > hardware::buffers::get_prng_buffer_size(get_device(), kernelParameters->getUseSameRndNumbers())) {
			*gs = hardware::buffers::get_prng_buffer_size(get_device(), kernelParameters->getUseSameRndNumbers());
		}
	}

	//Query specific sizes for kernels if needed
	if(kernel == scalar_product_eoprec || kernel == scalar_product || kernel == global_squarenorm || kernel == global_squarenorm_eoprec) {
		if(*ls > 64) {
			*ls = 64;
			*num_groups = (*gs)/(*ls);
		}
		return;
	}
}

void hardware::code::Spinors::convert_from_eoprec_device(const hardware::buffers::Spinor * in1, const hardware::buffers::Spinor * in2, const hardware::buffers::Plain<spinor> * out) const
{
	using namespace hardware::buffers;

	const size_4 mem_size = get_device()->getLocalLatticeMemoryExtents();

	// check buffer sizes
	const size_t in_size = in1->get_elements();
	if(in_size != get_eoprec_spinorfieldsize(mem_size) || in2->get_elements() != in_size) {
		throw std::invalid_argument("input buffers must be of size VOL4D / 2");
	}
	const size_t out_size = out->get_elements();
	if(out_size != get_spinorfieldsize(mem_size)) {
		throw std::invalid_argument("output buffer must be of size VOL4D");
	}

	const hardware::buffers::Buffer * tmp1, * tmp2;
	if(in1->is_soa()) {
		Plain<spinor> * tmp = new Plain<spinor>(in_size, get_device());
		convertSpinorfieldFromSOA_eo_device(tmp, in1);
		tmp1 = tmp;
	} else {
		tmp1 = in1;
	}
	if(in2->is_soa()) {
		Plain<spinor> * tmp = new Plain<spinor>(in_size, get_device());
		convertSpinorfieldFromSOA_eo_device(tmp, in2);
		tmp2 = tmp;
	} else {
		tmp2 = in2;
	}

	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(convert_from_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(convert_from_eoprec, 0, sizeof(cl_mem), tmp1->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(convert_from_eoprec, 1, sizeof(cl_mem), tmp2->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(convert_from_eoprec, 2, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(convert_from_eoprec, gs2, ls2);

	if(tmp1 != in1) {
		delete tmp1;
	}
	if(tmp2 != in2) {
		delete tmp2;
	}
}

void hardware::code::Spinors::convert_to_eoprec_device(const hardware::buffers::Spinor * out1, const hardware::buffers::Spinor * out2, const hardware::buffers::Plain<spinor> * in) const
{
	using namespace hardware::buffers;

	const size_4 mem_size = get_device()->getLocalLatticeMemoryExtents();

	// check buffer sizes
	const size_t out_size = out1->get_elements();
	if(out_size != get_eoprec_spinorfieldsize(mem_size) || out2->get_elements() != out_size) {
		throw std::invalid_argument("output buffers must be of size VOL4D / 2");
	}
	const size_t in_size = in->get_elements();
	if(in_size != get_spinorfieldsize(mem_size)) {
		throw std::invalid_argument("input buffer must be of size VOL4D");
	}

	const hardware::buffers::Buffer * tmp1, * tmp2;
	if(out1->is_soa()) {
		tmp1 = new Plain<spinor>(out_size, get_device());
	} else {
		tmp1 = out1;
	}
	if(out2->is_soa()) {
		tmp2 = new Plain<spinor>(out_size, get_device());
	} else {
		tmp2 = out2;
	}

	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(convert_to_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(convert_to_eoprec, 0, sizeof(cl_mem), tmp1->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(convert_to_eoprec, 1, sizeof(cl_mem), tmp2->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(convert_to_eoprec, 2, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(convert_to_eoprec , gs2, ls2);

	if(out1->is_soa()) {
		convertSpinorfieldToSOA_eo_device(out1, static_cast<const Plain<spinor> *>(tmp1));
		delete tmp1;
	}
	if(out2->is_soa()) {
		convertSpinorfieldToSOA_eo_device(out2, static_cast<const Plain<spinor> *>(tmp2));
		delete tmp2;
	}
}

//BLAS-functions
void hardware::code::Spinors::saxpy_device(const hardware::buffers::Plain<spinor> * x, const hardware::buffers::Plain<spinor> * y, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Plain<spinor> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxpy, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxpy, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy, 1, sizeof(cl_mem), y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy, 2, sizeof(cl_mem), alpha->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy, 3, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(saxpy , gs2, ls2);
}

void hardware::code::Spinors::saxpy_device(const hardware::buffers::Plain<spinor> * x, const hardware::buffers::Plain<spinor> * y, const hmc_complex alpha, const hardware::buffers::Plain<spinor> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxpy_arg, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxpy_arg, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_arg, 1, sizeof(cl_mem), y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_arg, 2, sizeof(hmc_float), &alpha.re);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_arg, 3, sizeof(hmc_float), &alpha.im);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_arg, 4, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(saxpy_arg , gs2, ls2);
}

void hardware::code::Spinors::sax_device(const hardware::buffers::Plain<spinor> * x, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Plain<spinor> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(sax, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(sax, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(sax, 1, sizeof(cl_mem), alpha->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(sax, 2, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(sax , gs2, ls2);
}

void hardware::code::Spinors::set_spinorfield_cold_device(const hardware::buffers::Plain<spinor> * inout) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(set_spinorfield_cold, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(set_spinorfield_cold, 0, sizeof(cl_mem), inout->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(set_spinorfield_cold , gs2, ls2);
}

void hardware::code::Spinors::set_eoprec_spinorfield_cold_device(const hardware::buffers::Spinor * inout) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(set_eoprec_spinorfield_cold, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(set_eoprec_spinorfield_cold, 0, sizeof(cl_mem), inout->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(set_eoprec_spinorfield_cold , gs2, ls2);
}

void hardware::code::Spinors::saxpy_eoprec_device(const hardware::buffers::Spinor * x, const hardware::buffers::Spinor * y, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Spinor * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxpy_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxpy_eoprec, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_eoprec, 1, sizeof(cl_mem), y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_eoprec, 2, sizeof(cl_mem), alpha->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_eoprec, 3, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( saxpy_eoprec, gs2, ls2);
}

void hardware::code::Spinors::saxpy_eoprec_device(const hardware::buffers::Spinor * x, const hardware::buffers::Spinor * y, const hmc_complex alpha, const hardware::buffers::Spinor * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxpy_arg_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxpy_arg_eoprec, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_arg_eoprec, 1, sizeof(cl_mem), y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_arg_eoprec, 2, sizeof(hmc_float), &alpha.re);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_arg_eoprec, 3, sizeof(hmc_float), &alpha.im);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_arg_eoprec, 4, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( saxpy_arg_eoprec, gs2, ls2);
}

void hardware::code::Spinors::sax_eoprec_device(const hardware::buffers::Spinor * x, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Spinor * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(sax_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(sax_eoprec, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(sax_eoprec, 1, sizeof(cl_mem), alpha->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(sax_eoprec, 2, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( sax_eoprec, gs2, ls2);
}

void hardware::code::Spinors::saxsbypz_device(const hardware::buffers::Plain<spinor> * x, const hardware::buffers::Plain<spinor> * y, const hardware::buffers::Plain<spinor> * z, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Plain<hmc_complex> * beta, const hardware::buffers::Plain<spinor> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxsbypz, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxsbypz, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxsbypz, 1, sizeof(cl_mem), y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxsbypz, 2, sizeof(cl_mem), z->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxsbypz, 3, sizeof(cl_mem), alpha->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxsbypz, 4, sizeof(cl_mem), beta->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxsbypz, 5, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(saxsbypz, gs2, ls2);
}

void hardware::code::Spinors::saxsbypz_eoprec_device(const hardware::buffers::Spinor * x, const hardware::buffers::Spinor * y, const hardware::buffers::Spinor * z, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Plain<hmc_complex> * beta, const hardware::buffers::Spinor * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxsbypz_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(saxsbypz_eoprec, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxsbypz_eoprec, 1, sizeof(cl_mem), y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxsbypz_eoprec, 2, sizeof(cl_mem), z->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxsbypz_eoprec, 3, sizeof(cl_mem), alpha->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxsbypz_eoprec, 4, sizeof(cl_mem), beta->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxsbypz_eoprec, 5, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(saxsbypz_eoprec, gs2, ls2);
}

void hardware::code::Spinors::set_complex_to_scalar_product_device(const hardware::buffers::Plain<spinor> * a, const hardware::buffers::Plain<spinor> * b, const hardware::buffers::Plain<hmc_complex> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(scalar_product, &ls2, &gs2, &num_groups);

	hardware::buffers::Plain<hmc_complex> tmp(num_groups, get_device());

	//set arguments
	int clerr = clSetKernelArg(scalar_product, 0, sizeof(cl_mem), a->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product, 1, sizeof(cl_mem), b->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product, 2, sizeof(cl_mem), tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product, 3, sizeof(hmc_complex) * ls2, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(scalar_product , gs2, ls2);


	clerr = clSetKernelArg(scalar_product_reduction, 0, sizeof(cl_mem), tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product_reduction, 1, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product_reduction, 2, sizeof(cl_uint), &num_groups);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( scalar_product_reduction, gs2, ls2);
}

void hardware::code::Spinors::set_complex_to_scalar_product_eoprec_device(const hardware::buffers::Spinor * a, const hardware::buffers::Spinor * b, const hardware::buffers::Plain<hmc_complex> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(scalar_product_eoprec, &ls2, &gs2, &num_groups);

	assert(scalar_product_buf->get_elements() == num_groups);

	//set arguments
	int clerr = clSetKernelArg(scalar_product_eoprec, 0, sizeof(cl_mem), a->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product_eoprec, 1, sizeof(cl_mem), b->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product_eoprec, 2, sizeof(cl_mem), scalar_product_buf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product_eoprec, 3, sizeof(hmc_complex) * ls2, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( scalar_product_eoprec, gs2, ls2);

	clerr = clSetKernelArg(scalar_product_reduction, 0, sizeof(cl_mem), scalar_product_buf->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product_reduction, 1, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(scalar_product_reduction, 2, sizeof(cl_uint), &num_groups);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(scalar_product_reduction, gs2, ls2);
}

void hardware::code::Spinors::global_squarenorm_reduction(const hardware::buffers::Plain<hmc_float> * out, const hardware::buffers::Plain<hmc_float> * tmp_buf) const
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

void hardware::code::Spinors::set_float_to_global_squarenorm_device(const hardware::buffers::Plain<spinor> * a, const hardware::buffers::Plain<hmc_float> * out) const
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

void hardware::code::Spinors::set_float_to_global_squarenorm_eoprec_device(const hardware::buffers::Spinor * a, const hardware::buffers::Plain<hmc_float> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(global_squarenorm_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	hardware::buffers::Plain<hmc_float> tmp(num_groups, get_device());

	int clerr = clSetKernelArg(global_squarenorm_eoprec, 0, sizeof(cl_mem), a->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	//CP: these do not have to be args of the function since they are global objects to the class opencl??
	clerr = clSetKernelArg(global_squarenorm_eoprec, 1, sizeof(cl_mem), tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(global_squarenorm_eoprec, 2, sizeof(hmc_float) * ls2, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( global_squarenorm_eoprec, gs2, ls2);

	global_squarenorm_reduction(out, &tmp);
}

void hardware::code::Spinors::set_zero_spinorfield_device(const hardware::buffers::Plain<spinor> * x) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(set_zero_spinorfield, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(set_zero_spinorfield, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(set_zero_spinorfield , gs2, ls2);
}

void hardware::code::Spinors::set_zero_spinorfield_eoprec_device(const hardware::buffers::Spinor * x) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(set_zero_spinorfield_eoprec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(set_zero_spinorfield_eoprec, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( set_zero_spinorfield_eoprec, gs2, ls2);
}

void hardware::code::Spinors::convertSpinorfieldToSOA_eo_device(const hardware::buffers::Spinor * out, const hardware::buffers::Plain<spinor> * in) const
{
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(convertSpinorfieldToSOA_eo, &ls2, &gs2, &num_groups);

	//set arguments
	int clerr = clSetKernelArg(convertSpinorfieldToSOA_eo, 0, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(convertSpinorfieldToSOA_eo, 1, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(convertSpinorfieldToSOA_eo, gs2, ls2);
}

void hardware::code::Spinors::convertSpinorfieldFromSOA_eo_device(const hardware::buffers::Plain<spinor> * out, const hardware::buffers::Spinor * in) const
{
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(convertSpinorfieldFromSOA_eo, &ls2, &gs2, &num_groups);

	//set arguments
	int clerr = clSetKernelArg(convertSpinorfieldFromSOA_eo, 0, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(convertSpinorfieldFromSOA_eo, 1, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(convertSpinorfieldFromSOA_eo, gs2, ls2);
}

// merged kernel calls
void hardware::code::Spinors::saxpy_AND_squarenorm_eo_device(const hardware::buffers::Spinor * x, const hardware::buffers::Spinor * y, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Spinor * out, const hardware::buffers::Plain<hmc_complex> * sq_out) const
{
	hmc_complex zero = {0.f, 0.f};
	sq_out->load(&zero);

	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(saxpy_AND_squarenorm_eo, &ls2, &gs2, &num_groups);

	//init local buffer for reduction
	hardware::buffers::Plain<hmc_float> tmp(num_groups, get_device());

	//set arguments
	int clerr = clSetKernelArg(saxpy_AND_squarenorm_eo, 0, sizeof(cl_mem), x->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_AND_squarenorm_eo, 1, sizeof(cl_mem), y->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_AND_squarenorm_eo, 2, sizeof(cl_mem), alpha->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_AND_squarenorm_eo, 3, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_AND_squarenorm_eo, 4, sizeof(cl_mem), tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(saxpy_AND_squarenorm_eo, 5, sizeof(hmc_float) * ls2, static_cast<void*>(nullptr));
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel( saxpy_AND_squarenorm_eo, gs2, ls2);

	clerr = clSetKernelArg(_global_squarenorm_reduction, 0, sizeof(cl_mem), sq_out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(_global_squarenorm_reduction, 1, sizeof(cl_mem), tmp);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	cl_uint elems = tmp.get_elements();
	clerr = clSetKernelArg(_global_squarenorm_reduction, 2, sizeof(cl_uint), &elems);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(_global_squarenorm_reduction, 1, 1);
}

size_t hardware::code::Spinors::get_read_write_size(const std::string& in) const
{
	//Depending on the compile-options, one has different sizes...
	size_t D = kernelParameters->getFloatSize();
	size_t S = kernelParameters->getSpinorFieldSize();
	size_t Seo = kernelParameters->getEoprecSpinorFieldSize();
	//factor for complex numbers
	int C = 2;
	//this is the same as in the function above
	//NOTE: 1 spinor has NC*NDIM = 12 complex entries
	if (in == "generate_gaussian_spinorfield") {
		//this kernel writes 1 spinor
		return ( 12 * C ) * D * S;
	}
	if (in == "generate_gaussian_spinorfield_eo") {
		//this kernel writes 1 spinor
		return ( 12 * C ) * D * Seo;
	}
	if (in == "set_spinorfield_cold") {
		//this kernel writes 1 spinor
		return C * 12 * D * S;
	}
	if (in == "set_eoprec_spinorfield_cold") {
		//this kernel writes 1 spinor
		return C * 12 * D * Seo;
	}
	if (in == "convert_from_eoprec") {
		//this kernel reads 2 spinor and writes 2 spinors per site
		///@todo is this right??
		return 2 * 2 * C * 12 * D * Seo;
	}
	if (in == "convert_to_eoprec") {
		//this kernel reads 2 spinor and writes 2 spinors per site
		return 2 * 2 * C * 12 * D * Seo;
	}
	if (in == "saxpy") {
		//this kernel reads 2 spinor, 2 complex number and writes 1 spinor per site
		return C * D * S * (12 * (2 + 1) + 2);
	}
	if (in == "sax") {
		//this kernel reads 1 spinor, 1 complex number and writes 1 spinor per site
		return C * D * S * (12 * (1 + 1) + 1);
	}
	if (in == "saxsbypz") {
		//this kernel reads 3 spinor, 2 complex number and writes 1 spinor per site
		return C * D * S * (12 * (3 + 1) + 2);
	}
	if (in == "set_zero_spinorfield") {
		//this kernel writes 1 spinor
		return C * 12 * D * S;
	}
	if (in == "saxpy_eoprec") {
		//this kernel reads 2 spinor, 1 complex number and writes 1 spinor per site
		return C * D * Seo * (12 * (2 + 1) + 1);
	}
	if (in == "sax_eoprec") {
		//this kernel reads 1 spinor, 1 complex number and writes 1 spinor per site
		return C * D * Seo * (12 * (1 + 1) + 1);
	}
	if (in == "saxsbypz_eoprec") {
		//this kernel reads 3 spinor, 2 complex number and writes 1 spinor per site
		return C * D * Seo * (12 * (3 + 1) + 2);
	}
	if (in == "set_zero_spinorfield_eoprec") {
		//this kernel writes 1 spinor
		return C * 12 * D * Seo;
	}
	if (in == "scalar_product") {
		//this kernel reads 2 spinors and writes 1 complex number
		/// @NOTE: here, the local reduction is not taken into account
		return C * D * S * ( 2 * 12  + 1 );
	}
	if (in == "scalar_product_reduction") {
		//this kernel reads NUM_GROUPS complex numbers and writes 1 complex number
		//query work-sizes for kernel to get num_groups
		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(scalar_product_reduction, &ls2, &gs2, &num_groups);
		return C * D * (num_groups + 1);
	}
	if (in == "global_squarenorm") {
		//this kernel reads 1 spinor and writes 1 real number
		/// @NOTE: here, the local reduction is not taken into account
		return D * S * (C * 12  + 1 );
	}
	if (in == "global_squarenorm_reduction") {
		//this kernel reads NUM_GROUPS real numbers and writes 1 real number
		//query work-sizes for kernel to get num_groups
		size_t ls2, gs2;
		cl_uint num_groups;
		this->get_work_sizes(scalar_product_reduction, &ls2, &gs2, &num_groups);
		return D * (num_groups + 1);
	}
	if (in == "scalar_product_eoprec") {
		//this kernel reads 2 spinors and writes 1 complex number
		/// @NOTE: here, the local reduction is not taken into account
		return C * D * Seo * ( 2 * 12  + 1 );
	}
	if (in == "global_squarenorm_eoprec") {
		//this kernel reads 1 spinor and writes 1 real number
		/// @NOTE: here, the local reduction is not taken into account
		return D * Seo * (C * 12  + 1 );
	}
	if(in == "convertSpinorfieldToSOA_eo") {
		return 2 * Seo * 24 * D;
	}
	if(in == "convertSpinorfieldFromSOA_eo") {
		return 2 * Seo * 24 * D;
	}
	//merged kernels
	if (in == "saxpy_AND_squarenorm_eo") {
		//the saxpy kernel reads 2 spinor, 2 complex number and writes 1 spinor per site
		//the squarenorm kernel reads 1 spinor and writes 1 real number
		// with the merging, the reading falls away
		/// @NOTE: here, the local reduction is not taken into account
		return C * D * Seo * (12 * (2 + 1) + 2)    +  D * Seo * ( 1 );
	}
	return 0;
}

uint64_t hardware::code::Spinors::get_flop_size(const std::string& in) const
{
	uint64_t S = kernelParameters->getSpinorFieldSize();
	uint64_t Seo = kernelParameters->getEoprecSpinorFieldSize();
	//this is the same as in the function above
	if (in == "generate_gaussian_spinorfield") {
		//this kernel performs 12 multiplications per site
		///@todo ? I did not count the gaussian normal pair production, which is very complicated...
		return 12 * S;
	}
	if (in == "generate_gaussian_spinorfield_eo") {
		//this kernel performs 12 multiplications per site
		///@todo ? I did not count the gaussian normal pair production, which is very complicated...
		return 12 * Seo;
	}
	if (in == "set_spinorfield_cold") {
		//this kernel performs 1. / sqrt((12.f * VOL4D)) and real_multiply_spinor for each site
		return S * ( 3 + 24);
	}
	if (in == "set_eoprec_spinorfield_cold") {
		//this kernel performs 1. / sqrt((12.f * VOL4D)) and real_multiply_spinor for each site
		return Seo * ( 3 + 24);
	}
	if (in == "convert_from_eoprec") {
		//this kernel does not perform any flop, he just copies memory
		return 0;
	}
	if (in == "convert_to_eoprec") {
		//this kernel does not perform any flop, he just copies memory
		return 0;
	}
	if (in == "saxpy") {
		//this kernel performs on each site spinor_times_complex and spinor_add
		return S * (NDIM * NC * ( getFlopComplexMult() + 2) );
	}
	if (in == "sax") {
		//this kernel performs on each site spinor_times_complex
		return S * (NDIM * NC * ( getFlopComplexMult() ) );
	}
	if (in == "saxsbypz") {
		//this kernel performs on each 2 * site spinor_times_complex and 2 * spinor_add
		return S * (NDIM * NC * 2 * ( getFlopComplexMult() + 2) );
	}
	if (in == "set_zero_spinorfield") {
		//this kernel does not do any flop
		return 0;
	}
	if (in == "saxpy_eoprec") {
		//this kernel performs on each site spinor_times_complex and spinor_add
		return Seo * (NDIM * NC * ( getFlopComplexMult() + 2) );
	}
	if (in == "sax_eoprec") {
		//this kernel performs on each site spinor_times_complex
		return Seo * (NDIM * NC * ( getFlopComplexMult() ) );
	}
	if (in == "saxsbypz_eoprec") {
		//this kernel performs on each 2 * site spinor_times_complex and 2 * spinor_add
		return Seo * (NDIM * NC * 2 * ( getFlopComplexMult() + 2) );
	}
	if (in == "set_zero_spinorfield_eoprec") {
		//this kernel does not do any flop
		return 0;
	}
	if (in == "scalar_product") {
		//this kernel performs spinor*spinor on each site and then adds S-1 complex numbers
		return S * getFlopSpinorTimesSpinor() + (S - 1) * 2;
	}
	if (in == "scalar_product_reduction") {
		return module_metric_not_implemented<uint64_t>();
	}
	if (in == "global_squarenorm") {
		//this kernel performs spinor_squarenorm on each site and then adds S-1 complex numbers
		return S * getFlopSpinorSquareNorm() + (S - 1) * 2;
	}
	if (in == "global_squarenorm_reduction") {
		return module_metric_not_implemented<uint64_t>();
	}
	if (in == "scalar_product_eoprec") {
		//this kernel performs spinor*spinor on each site and then adds S-1 complex numbers
		return Seo * getFlopSpinorTimesSpinor() + (Seo - 1) * 2;
	}
	if (in == "global_squarenorm_eoprec") {
		//this kernel performs spinor_squarenorm on each site and then adds S-1 complex numbers
		return Seo * getFlopSpinorSquareNorm() + (Seo - 1) * 2;
	}
	//merged kernels
	if (in == "saxpy_AND_squarenorm_eo") {
		//the saxpy kernel performs on each site spinor_times_complex and spinor_add
		//the squarenorm kernel performs spinor_squarenorm on each site and then adds S-1 complex numbers
		return Seo * (NDIM * NC * ( getFlopComplexMult() + 2) )  + Seo * getFlopSpinorSquareNorm() + (S - 1) * 2;
	}


	return 0;
}

void hardware::code::Spinors::print_profiling(const std::string& filename, int number) const
{
	Opencl_Module::print_profiling(filename, number);
	Opencl_Module::print_profiling(filename, set_spinorfield_cold);
	Opencl_Module::print_profiling(filename, set_eoprec_spinorfield_cold);
	Opencl_Module::print_profiling(filename, convert_from_eoprec);
	Opencl_Module::print_profiling(filename, convert_to_eoprec);
	Opencl_Module::print_profiling(filename, saxpy);
	Opencl_Module::print_profiling(filename, sax);
	Opencl_Module::print_profiling(filename, saxsbypz);
	Opencl_Module::print_profiling(filename, set_zero_spinorfield);
	Opencl_Module::print_profiling(filename, saxpy_eoprec);
	Opencl_Module::print_profiling(filename, sax_eoprec);
	Opencl_Module::print_profiling(filename, saxsbypz_eoprec);
	Opencl_Module::print_profiling(filename, set_zero_spinorfield_eoprec);
	Opencl_Module::print_profiling(filename, scalar_product);
	Opencl_Module::print_profiling(filename, scalar_product_reduction);
	Opencl_Module::print_profiling(filename, global_squarenorm);
	Opencl_Module::print_profiling(filename, _global_squarenorm_reduction);
	Opencl_Module::print_profiling(filename, scalar_product_eoprec);
	Opencl_Module::print_profiling(filename, global_squarenorm_eoprec);
	Opencl_Module::print_profiling(filename, convertSpinorfieldToSOA_eo);
	Opencl_Module::print_profiling(filename, convertSpinorfieldFromSOA_eo);
	Opencl_Module::print_profiling(filename, saxpy_AND_squarenorm_eo);
	Opencl_Module::print_profiling(filename, generate_gaussian_spinorfield);
	Opencl_Module::print_profiling(filename, generate_gaussian_spinorfield_eo);
}

void hardware::code::Spinors::copy_to_eoprec_spinorfield_buffer(const hardware::buffers::Spinor * buf, const spinor * const source) const
{
	using namespace hardware::buffers;

	Plain<spinor> tmp(buf->get_elements(), get_device());
	tmp.load(source);
	convertSpinorfieldToSOA_eo_device(buf, &tmp);
}

hardware::code::Spinors::Spinors(const hardware::code::OpenClKernelParametersInterface& kernelParameters, const hardware::Device * device)
	: Opencl_Module(kernelParameters, device), generate_gaussian_spinorfield(0), generate_gaussian_spinorfield_eo(0), saxpy_AND_squarenorm_eo(0)
{
	fill_kernels();

	if(kernelParameters.getUseEo() ) {
		size_t foo1, foo2;
		cl_uint groups;
		this->get_work_sizes(scalar_product_eoprec, &foo1, &foo2, &groups);
		scalar_product_buf = new hardware::buffers::Plain<hmc_complex>(groups, get_device());
	} else {
		scalar_product_buf = nullptr;
	}
}

hardware::code::Spinors::~Spinors()
{
	if(scalar_product_buf) {
		delete scalar_product_buf;
	}
	clear_kernels();
}

ClSourcePackage hardware::code::Spinors::get_sources() const noexcept
{
	return basic_fermion_code;
}

void hardware::code::Spinors::generate_gaussian_spinorfield_device(const hardware::buffers::Plain<spinor> * in, const hardware::buffers::PRNGBuffer * prng) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(generate_gaussian_spinorfield, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(generate_gaussian_spinorfield, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(generate_gaussian_spinorfield, 1, sizeof(cl_mem), prng->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(generate_gaussian_spinorfield  , gs2, ls2);

	if(logger.beDebug()) {
		hardware::buffers::Plain<hmc_float> force_tmp(1, get_device());
		hmc_float resid;
		get_device()->getSpinorCode()->set_float_to_global_squarenorm_device(in, &force_tmp);
		force_tmp.dump(&resid);
		logger.debug() <<  "\tinit gaussian spinorfield:\t" << resid;
		if(resid != resid) {
			throw Print_Error_Message("calculation of gaussian spinorfield gave nan! Aborting...", __FILE__, __LINE__);
		}
	}
}

void hardware::code::Spinors::generate_gaussian_spinorfield_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::PRNGBuffer * prng) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(generate_gaussian_spinorfield_eo, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(generate_gaussian_spinorfield_eo, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(generate_gaussian_spinorfield_eo, 1, sizeof(cl_mem), prng->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(generate_gaussian_spinorfield_eo, gs2, ls2);

	if(logger.beDebug()) {
		hardware::buffers::Plain<hmc_float> force_tmp(1, get_device());
		hmc_float resid;
		get_device()->getSpinorCode()->set_float_to_global_squarenorm_eoprec_device(in, &force_tmp);
		force_tmp.dump(&resid);
		logger.debug() <<  "\tinit gaussian spinorfield energy:\t" << resid;
		if(resid != resid) {
			throw Print_Error_Message("calculation of gaussian spinorfield gave nan! Aborting...", __FILE__, __LINE__);
		}
	}

}

size_t hardware::code::get_spinorfieldsize(const size_4& params)
{
	return get_vol4d(params);
}

size_t hardware::code::get_eoprec_spinorfieldsize(const size_4& params)
{
	return get_spinorfieldsize(params) / 2;
}
