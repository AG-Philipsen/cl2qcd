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
	cl_int err = clEnqueueWriteBuffer(*device, cl_buffer, CL_TRUE, 0, bytes, array, 0, nullptr, nullptr);
	if(err) {
		throw hardware::OpenclException(err, "clEnqueueWriteBuffer", __FILE__, __LINE__);
	}
}

void hardware::buffers::Buffer::dump(void * array) const
{
	cl_int err = clEnqueueReadBuffer(*device, cl_buffer, CL_TRUE, 0, bytes, array, 0, nullptr, nullptr);
	if(err) {
		throw hardware::OpenclException(err, "clEnqueueReadBuffer", __FILE__, __LINE__);
	}
}

hardware::SynchronizationEvent hardware::buffers::Buffer::dump_async(void * array) const
{
	cl_event event_cl;
	cl_int err = clEnqueueReadBuffer(*device, cl_buffer, CL_FALSE, 0, bytes, array, 0, nullptr, &event_cl);
	if(err) {
		throw hardware::OpenclException(err, "clEnqueueReadBuffer", __FILE__, __LINE__);
	}

	const hardware::SynchronizationEvent event(event_cl);
	err = clReleaseEvent(event_cl);
	if(err) {
		throw hardware::OpenclException(err, "clReleaseEvent", __FILE__, __LINE__);
	}
	return event;
}

const cl_mem* hardware::buffers::Buffer::get_cl_buffer() const noexcept
{
	return &cl_buffer;
}

hardware::Device * hardware::buffers::Buffer::get_device() const noexcept
{
	return device;
}

void hardware::buffers::Buffer::copyData(const Buffer* orig) const
{
	if(this->bytes != orig->bytes) {
		throw std::invalid_argument("The source and destination buffer must be of equal size!");
	} else {
		/*
		 * Now we have to play with the device a little.
		 * It seems on AMD hardware the buffer copy thing either pretty much sucks or I am using it wrong.
		 */
		const std::string dev_name = device->get_name();
		if(this->bytes == 16 && (dev_name == "Cypress" || dev_name == "Cayman")) {
			logger.debug() << "Using an OpenCL kernel to copy 16 bytes on " << dev_name << '.';
			device->get_buffer_code()->copy_16_bytes(this, orig);
		} else {
			logger.debug() << "Using default OpenCL buffer copy method for " << this->bytes << " bytes on " << dev_name << '.';
			int err = clEnqueueCopyBuffer(device->get_queue(), orig->cl_buffer, this->cl_buffer, 0, 0, this->bytes, 0, 0, 0);
			if(err) {
				throw hardware::OpenclException(err, "clEnqueueCopyBuffer", __FILE__, __LINE__);
			}
		}
	}
}

void hardware::buffers::Buffer::clear() const
{
#ifdef CL_VERSION_1_2
	if(sizeof(hmc_complex_zero) % bytes) {
		cl_char foo = 0;
		cl_int err = clEnqueueFillBuffer(*device, cl_buffer, &foo, sizeof(foo), 0, bytes, 0, nullptr, nullptr);
		if(err) {
			throw hardware::OpenclException(err, "clEnqueueFillBuffer", __FILE__, __LINE__);
		}
	} else {
		cl_int err = clEnqueueFillBuffer(*device, cl_buffer, &hmc_complex_zero, sizeof(hmc_complex_zero), 0, bytes, 0, nullptr, nullptr);
		if(err) {
			throw hardware::OpenclException(err, "clEnqueueFillBuffer", __FILE__, __LINE__);
		}
	}
#else
	device->get_buffer_code()->clear(this);
#endif
}
