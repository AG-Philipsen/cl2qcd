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

#include "heatbath.hpp"

#include "../../host_functionality/logger.hpp"
#include "../device.hpp"
#include "flopUtilities.hpp"
#include "gaugefield.hpp"
#include "prng.hpp"
#include "spinors.hpp"

using namespace std;

void hardware::code::Heatbath::fill_kernels()
{
	ClSourcePackage sources = get_basic_sources() << get_device()->getPrngCode()->get_sources() << "operations_geometry.cl" << "operations_complex.h" << "operations_matrix_su3.cl" << "operations_matrix.cl" << "operations_gaugefield.cl";

	logger.debug() << "Creating Heatbath kernels...";
	
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
	
	logger.debug() << "Clearing Heatbath kernels...";

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
			*gs = std::min(kernelParameters->getSpatialLatticeVolume() * kernelParameters->getNt() / 2, hardware::buffers::get_prng_buffer_size(get_device(), kernelParameters->getUseSameRndNumbers() ));
		} else {
			*gs = std::min(get_device()->get_num_compute_units(), hardware::buffers::get_prng_buffer_size(get_device(), kernelParameters->getUseSameRndNumbers() ));
		}
		*ls = get_device()->get_preferred_local_thread_num();
		*num_groups = *gs / *ls;
	}
	return;
}

size_t hardware::code::Heatbath::get_read_write_size(const std::string& in) const
{
	//Depending on the compile-options, one has different sizes...
	size_t D = kernelParameters->getFloatSize();
	size_t R = kernelParameters->getMatSize();
	//factor for complex numbers
	int C = 2;
	const size_t VOL4D = kernelParameters->getLatticeVolume();
	//this is the same as in the function above
	if ( (in == "heatbath_even" ) || (in == "heatbath_odd") || (in == "overrelax_even") || (in == "overrelax_odd")) {
		//this kernel reads ingredients for 1 staple plus 1 su3matrix and writes 1 su3-matrix
		return VOL4D / 2 * C * D * R * (6 * (NDIM - 1) + 1 + 1 );
	}
	return 0;
}

uint64_t hardware::code::Heatbath::get_flop_size(const std::string& in) const
{
	const size_t VOL4D = kernelParameters->getLatticeVolume();
	//this is the same as in the function above
	///@NOTE: I do not distinguish between su3 and 3x3 matrices. This is a difference if one use e.g. REC12, but here one wants to have the "netto" flops for comparability.
	if ( (in == "heatbath_even" ) || (in == "heatbath_odd") ) {
		//this kernel calculates 1 staple (= 4*ND-1 su3_su3 + 2_ND-1 su3_add) plus NC*(2*su3_su3 80 flops for the su2 update)
		return VOL4D / 2 * (4 * (NDIM - 1) * getFlopSu3MatrixTimesSu3Matrix() + 2 * (NDIM - 1) * 18 + NC * (2 * getFlopSu3MatrixTimesSu3Matrix() + 80));
	}
	if ( (in == "overrelax_even") || (in == "overrelax_odd")) {
		//this kernel calculates 1 staple (= 4*ND-1 su3_su3 + 2_ND-1 su3_add) plus NC*(2*su3_su3 58 flops for the su2 update)
		return VOL4D / 2 * (4 * (NDIM - 1) * getFlopSu3MatrixTimesSu3Matrix() + 2 * (NDIM - 1) * 18 + NC * (2 * getFlopSu3MatrixTimesSu3Matrix() + 58));
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

hardware::code::Heatbath::Heatbath(const hardware::code::OpenClKernelParametersInterface& kernelParameters, const hardware::Device * device)
	: Opencl_Module(kernelParameters, device)
{
	fill_kernels();
}

hardware::code::Heatbath::~Heatbath()
{
	clear_kernels();
}
