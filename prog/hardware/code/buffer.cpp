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

hardware::code::Buffer::~Buffer()
{
	clReleaseKernel(_copy_16_bytes);
}
