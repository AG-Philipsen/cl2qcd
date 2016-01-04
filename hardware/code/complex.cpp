/** @file
 * Implementation of the hardware::code::Complex class
 *
 * Copyright (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.phys.uni-frankfurt.de>
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

#include "complex.hpp"

#include "../../host_functionality/logger.hpp"
#include "../device.hpp"
#include "gaugefield.hpp"
#include <string>
#include "flopUtilities.hpp"

void hardware::code::Complex::fill_kernels()
{
	basic_complex_code = ClSourcePackage("-I " + std::string(SOURCEDIR) + " -D _INKERNEL_" + ((kernelParameters->getPrecision() == 64) ? (std::string(" -D _USEDOUBLEPREC_") + " -D _DEVICE_DOUBLE_EXTENSION_KHR_") : "")) << "types.h" << "operations_complex.h";
	
	logger.debug() << "Creating Complex kernels...";
	
	convert = createKernel("convert_float_to_complex") << basic_complex_code << "complex_convert.cl";
	ratio = createKernel("complex_ratio") << basic_complex_code << "complex_ratio.cl";
	product = createKernel("complex_product") << basic_complex_code << "complex_product.cl";
	sum = createKernel("complex_sum") << basic_complex_code << "complex_sum.cl";
	difference = createKernel("complex_subtraction") << basic_complex_code << "complex_subtraction.cl";
}

void hardware::code::Complex::clear_kernels()
{
	cl_int clerr = CL_SUCCESS;

	logger.debug() << "Clearing Complex kernels...";
	
	clerr = clReleaseKernel(convert);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(ratio);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(product);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(sum);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(difference);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
}


void hardware::code::Complex::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
	Opencl_Module::get_work_sizes(kernel, ls, gs, num_groups);

	//Query specific sizes for kernels if needed
	std::string kernelname = get_kernel_name(kernel);
	if(kernelname.compare("convert_float_to_complex") == 0) {
		*ls = 1;
		*gs = 1;
		*num_groups = 1;
	}
	if(kernelname.compare("complex_ratio") == 0) {
		*ls = 1;
		*gs = 1;
		*num_groups = 1;
	}
	if(kernelname.compare("complex_product") == 0) {
		*ls = 1;
		*gs = 1;
		*num_groups = 1;
	}
	if(kernelname.compare("complex_sum") == 0) {
		*ls = 1;
		*gs = 1;
		*num_groups = 1;
	}
	if(kernelname.compare("complex_subtraction") == 0) {
		*ls = 1;
		*gs = 1;
		*num_groups = 1;
	}
}

void hardware::code::Complex::set_complex_to_float_device(const hardware::buffers::Plain<hmc_float> * in, const hardware::buffers::Plain<hmc_complex> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(convert, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(convert, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(convert, 1, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(convert, gs2, ls2);
}

void hardware::code::Complex::set_complex_to_ratio_device(const hardware::buffers::Plain<hmc_complex> * a, const hardware::buffers::Plain<hmc_complex> * b, const hardware::buffers::Plain<hmc_complex> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(ratio, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(ratio, 0, sizeof(cl_mem), a->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(ratio, 1, sizeof(cl_mem), b->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(ratio, 2, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(ratio, gs2, ls2);
}

void hardware::code::Complex::set_complex_to_product_device(const hardware::buffers::Plain<hmc_complex> * a, const hardware::buffers::Plain<hmc_complex> * b, const hardware::buffers::Plain<hmc_complex> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(product, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(product, 0, sizeof(cl_mem), a->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(product, 1, sizeof(cl_mem), b->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(product, 2, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(product, gs2, ls2);
}

void hardware::code::Complex::set_complex_to_sum_device(const hardware::buffers::Plain<hmc_complex> * a, const hardware::buffers::Plain<hmc_complex> * b, const hardware::buffers::Plain<hmc_complex> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(sum, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(sum, 0, sizeof(cl_mem), a->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(sum, 1, sizeof(cl_mem), b->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(sum, 2, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(sum, gs2, ls2);
}

void hardware::code::Complex::set_complex_to_difference_device(const hardware::buffers::Plain<hmc_complex> * a, const hardware::buffers::Plain<hmc_complex> * b, const hardware::buffers::Plain<hmc_complex> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(difference, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(difference, 0, sizeof(cl_mem), a->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(difference, 1, sizeof(cl_mem), b->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(difference, 2, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(difference, gs2, ls2);
}

size_t hardware::code::Complex::get_read_write_size(const std::string& in) const
{
	size_t D = kernelParameters->getFloatSize();
	//factor for complex numbers
	int C = 2;
	if (in == "convert_float_to_complex") {
		//this kernel reads 1 float and writes 1 complex number
		return (C + 1) * D;
	}
	if (in == "complex_ratio" || in == "complex_product" || in == "complex_sum" || in == "complex_subtraction") {
		//this kernel reads 2 complex numbers and writes 1 complex number
		return C * D * (2 + 1);
	}

	logger.warn() << "No if entered in Complex::get_read_write_size, returning 0...";
	return 0;
}

uint64_t hardware::code::Complex::get_flop_size(const std::string& in) const
{
    if (in == "convert_float_to_complex") {
        return 0;
    }
	if (in == "complex_ratio") {
		return 11;
	}
	if (in == "complex_product") {
		return getFlopComplexMult();
	}
	if (in == "complex_sum") {
		return 2;
	}
	if (in == "complex_subtraction") {
		return 2;
	}
	
	logger.warn() << "No if entered in Complex::get_flop_size, returning 0...";
	return 0;
}

void hardware::code::Complex::print_profiling(const std::string& filename, int number) const
{
	Opencl_Module::print_profiling(filename, number);
	Opencl_Module::print_profiling(filename, ratio);
	Opencl_Module::print_profiling(filename, convert);
	Opencl_Module::print_profiling(filename, product);
	Opencl_Module::print_profiling(filename, sum);
	Opencl_Module::print_profiling(filename, difference);
}

hardware::code::Complex::Complex(const hardware::code::OpenClKernelParametersInterface& kernelParameters, const hardware::Device * device)
	: Opencl_Module(kernelParameters, device), convert(0), ratio(0), product(0), sum(0), difference(0)
{
	fill_kernels();
}

hardware::code::Complex::~Complex()
{
	clear_kernels();
}

ClSourcePackage hardware::code::Complex::get_sources() const noexcept
{
	return basic_complex_code;
}
