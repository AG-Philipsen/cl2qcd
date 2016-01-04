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

#include "kappa.hpp"

#include "../../host_functionality/logger.hpp"
#include "../device.hpp"
#include "gaugefield.hpp"

using namespace std;

void hardware::code::Kappa::fill_kernels()
{
	ClSourcePackage sources = get_basic_sources() << "operations_geometry.cl" << "operations_complex.h" << "operations_matrix_su3.cl" << "operations_matrix.cl" << "operations_gaugefield.cl";

	logger.debug() << "Creating TK clover kernels...";
	
	kappa_clover_gpu = createKernel("kappa_clover_gpu") << sources << "opencl_tk_kappa.cl";
}

void hardware::code::Kappa::run_kappa_clover(const hardware::buffers::Plain<hmc_float> * kappa, const hardware::buffers::SU3 * gaugefield, const hmc_float beta) const
{
	cl_int clerr = CL_SUCCESS;
	
	logger.debug() << "Clearing TK clover kernels...";

	size_t local_work_size;
	size_t global_work_size;
	cl_uint num_groups;

	get_work_sizes(kappa_clover_gpu, &local_work_size, &global_work_size, &num_groups);

	hardware::buffers::Plain<hmc_float> clmem_kappa_clover_buf_glob(num_groups, get_device());

	clerr = clSetKernelArg(kappa_clover_gpu, 0, sizeof(cl_mem), gaugefield->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(kappa_clover_gpu, 1, sizeof(hmc_float), &beta);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clKernelArg", __FILE__, __LINE__);

	clerr = clSetKernelArg(kappa_clover_gpu, 2, sizeof(cl_mem), kappa->get_cl_buffer());
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clKernelArg", __FILE__, __LINE__);

	get_device()->enqueue_kernel(kappa_clover_gpu, global_work_size, local_work_size);
}

void hardware::code::Kappa::clear_kernels()
{
	cl_int clerr = clReleaseKernel(kappa_clover_gpu);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clReleaseKernel", __FILE__, __LINE__);
}

hardware::code::Kappa::Kappa(const hardware::code::OpenClKernelParametersInterface& kernelParameters, const hardware::Device * device)
	: Opencl_Module(kernelParameters, device)
{
	fill_kernels();
}

hardware::code::Kappa::~Kappa()
{
	clear_kernels();
}
