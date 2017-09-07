/** @file
 * Implementation of the hardware::code::Real class
 *
 * Copyright (c) 2014 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
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

#include "real.hpp"

#include "../../host_functionality/logger.hpp"
#include "../device.hpp"

void hardware::code::Real::fill_kernels()
{
	basic_real_code = ClSourcePackage("-I " + std::string(SOURCEDIR) + " -D _INKERNEL_" + ((kernelParameters->getPrecision() == 64) ? (std::string(" -D _USEDOUBLEPREC_") + " -D _DEVICE_DOUBLE_EXTENSION_KHR_") : "")) << "types.h" << "operations_real.cl";
	
	logger.debug() << "Creating Real kernels...";
	
	//Setting operations kernel
	get_elem_vec = createKernel("get_elem_vector") << basic_real_code << "real_access_vector_element.cl";
	set_elem_vec = createKernel("set_elem_vector") << basic_real_code << "real_access_vector_element.cl";
	//Single operations kernels
	ratio = createKernel("real_ratio") << basic_real_code << "real_ratio.cl";
	product = createKernel("real_product") << basic_real_code << "real_product.cl";
	sum = createKernel("real_sum") << basic_real_code << "real_sum.cl";
	difference = createKernel("real_subtraction") << basic_real_code << "real_subtraction.cl";
	//Update cgm kernels
	update_alpha_cgm = createKernel("update_alpha_cgm") << basic_real_code << "update_alpha_cgm.cl";
	update_beta_cgm = createKernel("update_beta_cgm") << basic_real_code << "update_beta_cgm.cl";
	update_zeta_cgm = createKernel("update_zeta_cgm") << basic_real_code << "update_zeta_cgm.cl";
}

void hardware::code::Real::clear_kernels()
{
	cl_int clerr = CL_SUCCESS;

	logger.debug() << "Clearing Real kernels...";
	
	//Setting operations kernel
	clerr = clReleaseKernel(get_elem_vec);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(set_elem_vec);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	//Single operations kernels
	clerr = clReleaseKernel(ratio);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(product);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(sum);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(difference);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	//Update cgm kernels
	clerr = clReleaseKernel(update_alpha_cgm);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(update_beta_cgm);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
	clerr = clReleaseKernel(update_zeta_cgm);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
}


void hardware::code::Real::get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const
{
	Opencl_Module::get_work_sizes(kernel, ls, gs, num_groups);

	//Query specific sizes for kernels if needed
	std::string kernelname = get_kernel_name(kernel);
	if(kernelname.compare("get_elem_vector") == 0) {
		*ls = 1;
		*gs = 1;
		*num_groups = 1;
	}
	if(kernelname.compare("set_elem_vector") == 0) {
		*ls = 1;
		*gs = 1;
		*num_groups = 1;
	}
	if(kernelname.compare("real_ratio") == 0) {
		*ls = 1;
		*gs = 1;
		*num_groups = 1;
	}
	if(kernelname.compare("real_product") == 0) {
		*ls = 1;
		*gs = 1;
		*num_groups = 1;
	}
	if(kernelname.compare("real_sum") == 0) {
		*ls = 1;
		*gs = 1;
		*num_groups = 1;
	}
	if(kernelname.compare("real_subtraction") == 0) {
		*ls = 1;
		*gs = 1;
		*num_groups = 1;
	}
	if(kernelname.compare("update_alpha_cgm") == 0) {
		*ls = 16;
		*gs = 16;
		*num_groups = 1;
	}
	if(kernelname.compare("update_beta_cgm") == 0) {
		*ls = 16;
		*gs = 16;
		*num_groups = 1;
	}
	if(kernelname.compare("update_zeta_cgm") == 0) {
		*ls = 16;
		*gs = 16;
		*num_groups = 1;
	}
}

void hardware::code::Real::set_real_to_vector_element_device(const hardware::buffers::Plain<hmc_float> * in, const int index, const hardware::buffers::Plain<hmc_float> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(get_elem_vec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(get_elem_vec, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(get_elem_vec, 1, sizeof(cl_int), &index);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	clerr = clSetKernelArg(get_elem_vec, 2, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(get_elem_vec, gs2, ls2);
}

void hardware::code::Real::set_vector_element_to_real_device(const hardware::buffers::Plain<hmc_float> * in, const int index, const hardware::buffers::Plain<hmc_float> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(set_elem_vec, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(set_elem_vec, 0, sizeof(cl_mem), in->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(set_elem_vec, 1, sizeof(cl_int), &index);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	clerr = clSetKernelArg(set_elem_vec, 2, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(set_elem_vec, gs2, ls2);
}

void hardware::code::Real::set_real_to_ratio_device(const hardware::buffers::Plain<hmc_float> * a, const hardware::buffers::Plain<hmc_float> * b, const hardware::buffers::Plain<hmc_float> * out) const
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

void hardware::code::Real::set_real_to_product_device(const hardware::buffers::Plain<hmc_float> * a, const hardware::buffers::Plain<hmc_float> * b, const hardware::buffers::Plain<hmc_float> * out) const
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

void hardware::code::Real::set_real_to_sum_device(const hardware::buffers::Plain<hmc_float> * a, const hardware::buffers::Plain<hmc_float> * b, const hardware::buffers::Plain<hmc_float> * out) const
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

void hardware::code::Real::set_real_to_difference_device(const hardware::buffers::Plain<hmc_float> * a, const hardware::buffers::Plain<hmc_float> * b, const hardware::buffers::Plain<hmc_float> * out) const
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



void hardware::code::Real::update_zeta_cgm_device(const hardware::buffers::Plain<hmc_float> * zeta_prev, const hardware::buffers::Plain<hmc_float> * zeta_prev_prev, const hardware::buffers::Plain<hmc_float> * sbeta_prev, const hardware::buffers::Plain<hmc_float> * sbeta_pres, const hardware::buffers::Plain<hmc_float> * salpha_prev, const hardware::buffers::Plain<hmc_float> * sigma, const int numeq, const hardware::buffers::Plain<hmc_float> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(update_zeta_cgm, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(update_zeta_cgm, 0, sizeof(cl_mem), zeta_prev->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(update_zeta_cgm, 1, sizeof(cl_mem), zeta_prev_prev->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	clerr = clSetKernelArg(update_zeta_cgm, 2, sizeof(cl_mem), sbeta_prev->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	clerr = clSetKernelArg(update_zeta_cgm, 3, sizeof(cl_mem), sbeta_pres->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	clerr = clSetKernelArg(update_zeta_cgm, 4, sizeof(cl_mem), salpha_prev->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	clerr = clSetKernelArg(update_zeta_cgm, 5, sizeof(cl_mem), sigma->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	clerr = clSetKernelArg(update_zeta_cgm, 6, sizeof(int), &numeq);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	clerr = clSetKernelArg(update_zeta_cgm, 7, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(update_zeta_cgm, gs2, ls2);
}


void hardware::code::Real::update_beta_cgm_device(const hardware::buffers::Plain<hmc_float> * sbeta_pres, const hardware::buffers::Plain<hmc_float> * zeta_pres, const hardware::buffers::Plain<hmc_float> * zeta_prev, const int numeq, const hardware::buffers::Plain<hmc_float> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(update_beta_cgm, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(update_beta_cgm, 0, sizeof(cl_mem), sbeta_pres->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(update_beta_cgm, 1, sizeof(cl_mem), zeta_pres->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(update_beta_cgm, 2, sizeof(cl_mem), zeta_prev->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	clerr = clSetKernelArg(update_beta_cgm, 3, sizeof(int), &numeq);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	clerr = clSetKernelArg(update_beta_cgm, 4, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(update_beta_cgm, gs2, ls2);
}


void hardware::code::Real::update_alpha_cgm_device(const hardware::buffers::Plain<hmc_float> * salpha_pres, const hardware::buffers::Plain<hmc_float> * zeta_pres, const hardware::buffers::Plain<hmc_float> * beta_pres, const hardware::buffers::Plain<hmc_float> * zeta_prev, const hardware::buffers::Plain<hmc_float> * sbeta_pres, const int numeq, const hardware::buffers::Plain<hmc_float> * out) const
{
	//query work-sizes for kernel
	size_t ls2, gs2;
	cl_uint num_groups;
	this->get_work_sizes(update_alpha_cgm, &ls2, &gs2, &num_groups);
	//set arguments
	int clerr = clSetKernelArg(update_alpha_cgm, 0, sizeof(cl_mem), salpha_pres->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(update_alpha_cgm, 1, sizeof(cl_mem), zeta_pres->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	clerr = clSetKernelArg(update_alpha_cgm, 2, sizeof(cl_mem), beta_pres->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(update_alpha_cgm, 3, sizeof(cl_mem), zeta_prev->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(update_alpha_cgm, 4, sizeof(cl_mem), sbeta_pres->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(update_alpha_cgm, 5, sizeof(int), &numeq);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	
	clerr = clSetKernelArg(update_alpha_cgm, 6, sizeof(cl_mem), out->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(update_alpha_cgm, gs2, ls2);
}

size_t hardware::code::Real::get_read_write_size(const std::string& in) const
{
	size_t D = kernelParameters->getFloatSize();
	if (in == "ratio" || in == "product" || in == "sum" || in == "subtraction") {
		//this kernel reads 2 real numbers and writes 1 real number
		return D * (2 + 1);
	}
	if (in == "get_elem_vector" || in == "set_elem_vector") {
		//this kernel reads 1 real number and writes 1 real number
		return D * (1 + 1);
	}
	
	logger.warn() << "No if entered in Real::get_read_write_size, returning 0...";
	return 0;
}
uint64_t hardware::code::Real::get_flop_size(const std::string& in) const
{
	if (in == "ratio" || in == "product" || in == "sum" || in == "subtraction") {
		return 1;
	}
	if (in == "get_elem_vector" || in == "set_elem_vector") {
		return 0;
	}
	
	logger.warn() << "No if entered in Real::get_flop_size, returning 0...";
	return 0;
}

size_t hardware::code::Real::get_read_write_size_update(const std::string& in, const int numeq) const
{
	size_t D = kernelParameters->getFloatSize();
	if (in == "update_alpha_cgm") {
		if(numeq==0) 
			throw Print_Error_Message("get_read_write_size_update with update_alpha_cgm and numeq=0", __FILE__, __LINE__);
		//this kernel reads 3 vectors of float, 3 float, 1 int and writes 1 vector of float
		return (4 * numeq + 3 ) * D + sizeof(int);
	}
	if (in == "update_beta_cgm") {
		if(numeq==0) 
			throw Print_Error_Message("get_read_write_size_update with update_beta_cgm and numeq=0", __FILE__, __LINE__);
		//this kernel reads 2 vectors of float, 1 float, 1 int and writes 1 vector of float
		return (3 * numeq + 1 ) * D + sizeof(int);
	}
	if (in == "update_zeta_cgm") {
		if(numeq==0) 
			throw Print_Error_Message("get_read_write_size_update with update_zeta_cgm and numeq=0", __FILE__, __LINE__);
		//this kernel reads 3 vectors of float, 2 float, 1 int and writes 1 vector of float
		return (4 * numeq + 2 ) * D + sizeof(int);
	}
	return 0;
}

uint64_t hardware::code::Real::get_flop_size_update(const std::string& in, const int numeq) const
{
	if (in == "update_alpha_cgm") {
		if(numeq==0) 
			throw Print_Error_Message("get_flop_size_update with update_alpha_cgm and numeq=0", __FILE__, __LINE__);
		return 4 * numeq;
	}
	if (in == "update_beta_cgm") {
		if(numeq==0) 
			throw Print_Error_Message("get_flop_size_update with update_beta_cgm and numeq=0", __FILE__, __LINE__);
		return 2 * numeq;
	}
	if (in == "update_zeta_cgm") {
		if(numeq==0) 
			throw Print_Error_Message("get_flop_size_update with update_zeta_cgm and numeq=0", __FILE__, __LINE__);
		return 11 * numeq;
	}
	return 0;
}

void hardware::code::Real::print_profiling(const std::string& filename, int number) const
{
	Opencl_Module::print_profiling(filename, number);
	Opencl_Module::print_profiling(filename, get_elem_vec);
	Opencl_Module::print_profiling(filename, set_elem_vec);
	Opencl_Module::print_profiling(filename, ratio);
	Opencl_Module::print_profiling(filename, product);
	Opencl_Module::print_profiling(filename, sum);
	Opencl_Module::print_profiling(filename, difference);
	Opencl_Module::print_profiling(filename, update_alpha_cgm);
	Opencl_Module::print_profiling(filename, update_beta_cgm);
	Opencl_Module::print_profiling(filename, update_zeta_cgm);
}

hardware::code::Real::Real(const hardware::code::OpenClKernelParametersInterface& kernelParameters, const hardware::Device * device)
	: Opencl_Module(kernelParameters, device), get_elem_vec(0), set_elem_vec(0), ratio(0), product(0), sum(0),
	                                 difference(0), update_alpha_cgm(0), update_beta_cgm(0), update_zeta_cgm(0)
{
	fill_kernels();
}

hardware::code::Real::~Real()
{
	clear_kernels();
}

ClSourcePackage hardware::code::Real::get_sources() const noexcept
{
	return basic_real_code;
}
