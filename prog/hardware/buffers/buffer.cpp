/** @file
 * Implementation of the hardware::buffers::Buffer class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "buffer.hpp"

#include "../system.hpp"
#include "../../logger.hpp"

static cl_mem allocateBuffer(size_t bytes, cl_context context);

hardware::buffers::Buffer::Buffer(size_t bytes, hardware::Device * device)
	: bytes(bytes), cl_buffer(allocateBuffer(bytes, device->context))
{
	// nothing to do here, initialization complete
}

hardware::buffers::Buffer::~Buffer()
{
	clReleaseMemObject(cl_buffer);
}

static cl_mem allocateBuffer(size_t bytes, cl_context context)
{
	cl_int err;
	cl_mem cl_buffer = clCreateBuffer(context, 0, bytes, 0, &err);
	if(err) {
		throw hardware::OpenclException(err, "clCreateBuffer", __FILE__, __LINE__);
	}
	return cl_buffer;
}

hardware::buffers::Buffer::operator const cl_mem*() const
{
	return &cl_buffer;
}

size_t hardware::buffers::Buffer::get_bytes() const
{
	return bytes;
}
