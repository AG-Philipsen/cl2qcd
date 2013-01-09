/* @file
 * Implementation of the hardware::code::Buffer class
 */

#include "buffer.hpp"

#include <cassert>
#include "../device.hpp"

hardware::code::Buffer::Buffer(const meta::Inputparameters& params, hardware::Device * device)
	: Opencl_Module(params, device)
{
	_copy_16_bytes = createKernel("copy_16_bytes") << "buffer.cl";

	auto base_code = get_device()->get_gaugefield_code()->get_sources();
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
		get_device()->enqueue_kernel(_clear_bytes);
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
		get_device()->enqueue_kernel(_clear_float4);
	}
}

hardware::code::Buffer::~Buffer()
{
	clReleaseKernel(_copy_16_bytes);
}
