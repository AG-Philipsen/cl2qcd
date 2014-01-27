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
#include "../../meta/util.hpp"
#include "../device.hpp"
// #include "gaugefield.hpp"
// #include <string>

void hardware::code::Real::fill_kernels()
{
	basic_real_code = ClSourcePackage("-I " + std::string(SOURCEDIR) + " -D _INKERNEL_" + ((get_parameters().get_precision() == 64) ? (std::string(" -D _USEDOUBLEPREC_") + " -D _DEVICE_DOUBLE_EXTENSION_KHR_") : "")) << "types.h" << "operations_real.cl" << "real_algebra.cl";
	
	logger.debug() << "Creating Real kernels...";
	
	update_alpha_cgm = createKernel("update_alpha_cgm") << basic_real_code;
	update_beta_cgm = createKernel("update_beta_cgm") << basic_real_code;
	update_zeta_cgm = createKernel("update_zeta_cgm") << basic_real_code;
}

void hardware::code::Real::clear_kernels()
{
	cl_int clerr = CL_SUCCESS;

	logger.debug() << "Clearing Real kernels...";
	
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
	throw Print_Error_Message("get_read_write_size without methods!", __FILE__, __LINE__);
	return 0;
}
uint64_t hardware::code::Real::get_flop_size(const std::string& in) const
{
	throw Print_Error_Message("get_flop_size without methods!", __FILE__, __LINE__);
	return 0;
}

size_t hardware::code::Real::get_read_write_size_update(const std::string& in, const int numeq) const
{
	size_t D = meta::get_float_size(get_parameters());
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
	Opencl_Module::print_profiling(filename, update_alpha_cgm);
	Opencl_Module::print_profiling(filename, update_beta_cgm);
	Opencl_Module::print_profiling(filename, update_zeta_cgm);
}

hardware::code::Real::Real(const meta::Inputparameters& params, hardware::Device * device)
	: Opencl_Module(params, device), update_alpha_cgm(0), update_beta_cgm(0), update_zeta_cgm(0)
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
