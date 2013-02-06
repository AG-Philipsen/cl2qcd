/** @file
 * Implementation of the hardware::buffers::Buffer class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "buffer.hpp"

#include "../system.hpp"
#include "../../logger.hpp"
#include "../../crypto/md5.h"

static cl_mem allocateBuffer(size_t bytes, cl_context context, bool place_on_host);

hardware::buffers::Buffer::Buffer(size_t bytes, hardware::Device * device, bool place_on_host)
	: bytes(bytes), cl_buffer(allocateBuffer(bytes, device->context, place_on_host)), device(device)
{
	// nothing to do here, initialization complete
}

hardware::buffers::Buffer::~Buffer()
{
	clReleaseMemObject(cl_buffer);
}

static cl_mem allocateBuffer(size_t bytes, cl_context context, const bool place_on_host)
{
	cl_int err;
	const cl_mem_flags mem_flags = place_on_host ? CL_MEM_ALLOC_HOST_PTR : 0;
	cl_mem cl_buffer = clCreateBuffer(context, mem_flags, bytes, 0, &err);
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
			int err = clEnqueueCopyBuffer(device->get_queue(), orig->cl_buffer, this->cl_buffer, 0, 0, this->bytes, 0, nullptr, nullptr);
			if(err) {
				throw hardware::OpenclException(err, "clEnqueueCopyBuffer", __FILE__, __LINE__);
			}
		}
	}
}
void hardware::buffers::Buffer::copyDataBlock(const Buffer* orig, const size_t dest_offset, const size_t src_offset, const size_t bytes) const
{
	logger.debug() << "Copying " << bytes << " bytes from offset " << src_offset << " to offset " << dest_offset;
	logger.debug() << "Source buffer size: " << orig->bytes;
	logger.debug() << "Dest buffer size:   " << this->bytes;
	if(this->bytes < dest_offset + bytes || orig->bytes < src_offset + bytes) {
		throw std::invalid_argument("Copy range exceeds buffer size!");
	} else {
		int err = clEnqueueCopyBuffer(device->get_queue(), orig->cl_buffer, this->cl_buffer, src_offset, dest_offset, bytes, 0, nullptr, nullptr);
		if(err) {
			throw hardware::OpenclException(err, "clEnqueueCopyBuffer", __FILE__, __LINE__);
		}
	}
}

void hardware::buffers::Buffer::clear() const
{
#ifdef CL_VERSION_1_2
	device->synchronize();
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

std::string hardware::buffers::md5(const Buffer* buf)
{
	md5_t md5_state;
	md5_init(&md5_state);

	char* data = new char[buf->bytes];
	buf->dump(data);

	md5_process(&md5_state, data, buf->bytes);

	delete[] data;

	char sig[MD5_SIZE];
	md5_finish(&md5_state, sig);

	char res[33];
	md5_sig_to_string(sig, res, 33);

	return std::string(res);
}
