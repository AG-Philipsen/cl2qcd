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
	: bytes(bytes), cl_buffer(allocateBuffer(bytes, device->context)), device(device)
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

hardware::buffers::Buffer::operator const cl_mem*() const noexcept
{
	return &cl_buffer;
}

size_t hardware::buffers::Buffer::get_bytes() const noexcept
{
	return bytes;
}

void hardware::buffers::Buffer::load(const void * array) const
{
	cl_int err = clEnqueueWriteBuffer(*device, cl_buffer, CL_TRUE, 0, bytes, array, 0, 0, 0);
	if(err) {
		throw hardware::OpenclException(err, "clEnqueueWriteBuffer", __FILE__, __LINE__);
	}
}

void hardware::buffers::Buffer::dump(void * array) const
{
	cl_int err = clEnqueueReadBuffer(*device, cl_buffer, CL_TRUE, 0, bytes, array, 0, 0, 0);
	if(err) {
		throw hardware::OpenclException(err, "clEnqueueReadBuffer", __FILE__, __LINE__);
	}
}

const cl_mem* hardware::buffers::Buffer::get_cl_buffer() const noexcept
{
	return &cl_buffer;
}

hardware::Device * hardware::buffers::Buffer::get_device() const noexcept
{
	return device;
}
