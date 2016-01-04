/* @file
 * Implementation of the hardware::code::Buffer class
 *
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

#include "buffer.hpp"

#include <cassert>
#include "../device.hpp"

hardware::code::Buffer::Buffer(const hardware::code::OpenClKernelParametersInterface& kernelParameters, const hardware::Device * device)
	: Opencl_Module(kernelParameters, device)
{
	_copy_16_bytes = createKernel("copy_16_bytes") << "buffer.cl";

	auto base_code = get_device()->getGaugefieldCode()->get_sources();
	_clear_bytes = createKernel("clear_bytes") << base_code << "buffer.cl";
	_clear_float4 = createKernel("clear_float4") << base_code << "buffer.cl";
}

void hardware::code::Buffer::copy_16_bytes(const hardware::buffers::Buffer * dest, const hardware::buffers::Buffer * orig) const
{
	cl_int err = clSetKernelArg(_copy_16_bytes, 0, sizeof(cl_mem), dest->get_cl_buffer());
	if(err) {
		throw Opencl_Error(err, "clSetKernelArg", __FILE__, __LINE__);
	}

	err = clSetKernelArg(_copy_16_bytes, 1, sizeof(cl_mem), orig->get_cl_buffer());
	if(err) {
		throw Opencl_Error(err, "clSetKernelArg", __FILE__, __LINE__);
	}
	get_device()->enqueue_kernel(_copy_16_bytes, 1, 1);
}

void hardware::code::Buffer::clear(const hardware::buffers::Buffer * dest) const
{
	const cl_ulong bytes = dest->get_bytes();
	if(bytes % sizeof(cl_float4)) {
		cl_int err = clSetKernelArg(_clear_bytes, 0, sizeof(cl_mem), dest->get_cl_buffer());
		if(err) {
			throw Opencl_Error(err, "clSetKernelArg", __FILE__, __LINE__);
		}
		err = clSetKernelArg(_clear_bytes, 1, sizeof(cl_ulong), &bytes);
		if(err) {
			throw Opencl_Error(err, "clSetKernelArg", __FILE__, __LINE__);
		}
		get_device()->enqueueKernel(_clear_bytes);
	} else {
		cl_long elems = bytes / sizeof(cl_float4);
		cl_int err = clSetKernelArg(_clear_float4, 0, sizeof(cl_mem), dest->get_cl_buffer());
		if(err) {
			throw Opencl_Error(err, "clSetKernelArg", __FILE__, __LINE__);
		}
		err = clSetKernelArg(_clear_float4, 1, sizeof(cl_ulong), &elems);
		if(err) {
			throw Opencl_Error(err, "clSetKernelArg", __FILE__, __LINE__);
		}
		get_device()->enqueueKernel(_clear_float4);
	}
}

hardware::code::Buffer::~Buffer()
{
	clReleaseKernel(_copy_16_bytes);
}
