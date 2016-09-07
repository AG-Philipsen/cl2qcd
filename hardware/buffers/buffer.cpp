/** @file
 * Implementation of the hardware::buffers::Buffer class
 *
 * Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
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

#include "../system.hpp"
#include "../device.hpp"
#include "../../host_functionality/logger.hpp"
#include "../../crypto/md5.hpp"
#include "../code/buffer.hpp"

static cl_mem allocateBuffer(size_t bytes, cl_context context, bool place_on_host, cl_mem_flags extra_flags);
void memObjectReleased(cl_mem, void * user_data);
struct MemObjectAllocationTracer {
	size_t bytes;
	bool host;
	const hardware::Device * device;

	MemObjectAllocationTracer(size_t bytes, bool host, const hardware::Device * device)
	 : bytes(bytes), host(host), device(device) {
		device->markMemAllocated(host, bytes);
	 };

	~MemObjectAllocationTracer() {
		device->markMemReleased(host, bytes);
	}
};


hardware::buffers::Buffer::Buffer(const size_t bytes, const hardware::Device * device, const bool place_on_host, const cl_mem_flags extra_flags)
	: bytes(bytes), cl_buffer(allocateBuffer(bytes, device->context, place_on_host, extra_flags)), device(device)
{
	// notify device about allocation
	cl_int err = clSetMemObjectDestructorCallback(cl_buffer, memObjectReleased, new MemObjectAllocationTracer(bytes, place_on_host, device));
	if(err) {
		throw hardware::OpenclException(err, "clSetMemObjectDestructorCallback", __FILE__, __LINE__);
	}
}

hardware::buffers::Buffer::~Buffer()
{
	clReleaseMemObject(cl_buffer);
}

static cl_mem allocateBuffer(const size_t bytes, const cl_context context, const bool place_on_host, const cl_mem_flags extra_flags)
{
	cl_int err;
	const cl_mem_flags mem_flags = (place_on_host ? CL_MEM_ALLOC_HOST_PTR : 0) | extra_flags;
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

void hardware::buffers::Buffer::load(const void * array, size_t bytes, const size_t offset) const
{
	if(bytes == 0) {
		bytes = this->bytes;
	} else {
		if(bytes + offset > this->bytes) {
			logger.error() << "Writing outside buffer. Bytes: " << bytes << " Offset: " << offset << " Buffer size: " << this->bytes;
			throw std::out_of_range("You are loading to memory outside of the buffer.");
		}
	}
	logger.trace() << "clEnqueueWriteBuffer(...,...,...," << offset << ',' << bytes << ",...,0,nullptr,nullptr)";
	cl_int err = clEnqueueWriteBuffer(*device, cl_buffer, CL_TRUE, offset, bytes, array, 0, nullptr, nullptr);
	if(err) {
		throw hardware::OpenclException(err, "clEnqueueWriteBuffer", __FILE__, __LINE__);
	}
}

void hardware::buffers::Buffer::dump(void * array, size_t bytes, const size_t offset) const
{
	if(bytes == 0) {
		bytes = this->bytes;
	} else {
		if(bytes + offset > this->bytes) {
			logger.error() << "Reading outside buffer. Bytes: " << bytes << " Offset: " << offset << " Buffer size: " << this->bytes;
			throw std::out_of_range("You are reading from outside of the buffer.");
		}
	}
	logger.trace() << "clEnqueueReadBuffer(...,...,...," << offset << ',' << bytes << ",...,0,nullptr,nullptr)";
	cl_int err = clEnqueueReadBuffer(*device, cl_buffer, CL_TRUE, offset, bytes, array, 0, nullptr, nullptr);
	if(err) {
		throw hardware::OpenclException(err, "clEnqueueReadBuffer", __FILE__, __LINE__);
	}
}

void hardware::buffers::Buffer::load_rect(const void* src, const size_t *buffer_origin, const size_t *host_origin, const size_t *region, size_t buffer_row_pitch, size_t buffer_slice_pitch, size_t host_row_pitch, size_t host_slice_pitch) const
{
	cl_int err = clEnqueueWriteBufferRect(*device, cl_buffer, CL_TRUE, buffer_origin, host_origin, region, buffer_row_pitch, buffer_slice_pitch, host_row_pitch, host_slice_pitch, src, 0, nullptr, nullptr);
	if(err) {
		throw hardware::OpenclException(err, "clEnqueueReadBufferRect", __FILE__, __LINE__);
	}
}

void hardware::buffers::Buffer::dump_rect(void* src, const size_t *buffer_origin, const size_t *host_origin, const size_t *region, size_t buffer_row_pitch, size_t buffer_slice_pitch, size_t host_row_pitch, size_t host_slice_pitch) const
{
	cl_int err = clEnqueueReadBufferRect(*device, cl_buffer, CL_TRUE, buffer_origin, host_origin, region, buffer_row_pitch, buffer_slice_pitch, host_row_pitch, host_slice_pitch, src, 0, nullptr, nullptr);
	if(err) {
		throw hardware::OpenclException(err, "clEnqueueReadBufferRect", __FILE__, __LINE__);
	}
}

hardware::SynchronizationEvent hardware::buffers::Buffer::load_rectAsync(const void* src, const size_t *buffer_origin, const size_t *host_origin, const size_t *region, size_t buffer_row_pitch, size_t buffer_slice_pitch, size_t host_row_pitch, size_t host_slice_pitch, const hardware::SynchronizationEvent& event) const
{
	const cl_event * wait_event = nullptr;
	cl_uint num_wait_events = 0;
	if(event.raw()) {
		wait_event = &event.raw();
		num_wait_events = 1;
	}

	cl_event event_cl;
	cl_int err = clEnqueueWriteBufferRect(*device, cl_buffer, CL_FALSE, buffer_origin, host_origin, region, buffer_row_pitch, buffer_slice_pitch, host_row_pitch, host_slice_pitch, src, num_wait_events, wait_event, &event_cl);
	if(err) {
		throw hardware::OpenclException(err, "clEnqueueReadBufferRect", __FILE__, __LINE__);
	}

	const hardware::SynchronizationEvent new_event(event_cl);
	err = clReleaseEvent(event_cl);
	if(err) {
		throw hardware::OpenclException(err, "clReleaseEvent", __FILE__, __LINE__);
	}
	return new_event;
}

hardware::SynchronizationEvent hardware::buffers::Buffer::dump_rectAsync(void* dest, const size_t *buffer_origin, const size_t *host_origin, const size_t *region, size_t buffer_row_pitch, size_t buffer_slice_pitch, size_t host_row_pitch, size_t host_slice_pitch, const hardware::SynchronizationEvent& event) const
{
	const cl_event * wait_event = nullptr;
	cl_uint num_wait_events = 0;
	if(event.raw()) {
		wait_event = &event.raw();
		num_wait_events = 1;
	}

	cl_event event_cl;
	cl_int err = clEnqueueReadBufferRect(*device, cl_buffer, CL_FALSE, buffer_origin, host_origin, region, buffer_row_pitch, buffer_slice_pitch, host_row_pitch, host_slice_pitch, dest, num_wait_events, wait_event, &event_cl);
	if(err) {
		throw hardware::OpenclException(err, "clEnqueueReadBufferRect", __FILE__, __LINE__);
	}

	const hardware::SynchronizationEvent new_event(event_cl);
	err = clReleaseEvent(event_cl);
	if(err) {
		throw hardware::OpenclException(err, "clReleaseEvent", __FILE__, __LINE__);
	}
	return new_event;
}

hardware::SynchronizationEvent hardware::buffers::Buffer::load_async(const void * array) const
{
	cl_event event_cl;
	cl_int err = clEnqueueWriteBuffer(*device, cl_buffer, CL_FALSE, 0, bytes, array, 0, nullptr, &event_cl);
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

const hardware::Device * hardware::buffers::Buffer::get_device() const noexcept
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
			device->getBufferCode()->copy_16_bytes(this, orig);
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
	device->getBufferCode()->clear(this);
#endif
}

std::string hardware::buffers::md5(const Buffer* buf)
{
	std::vector<char> data(buf->bytes);
	buf->dump(data.data());
	return crypto::md5(std::string{begin(data), end(data)});
}

hardware::SynchronizationEvent hardware::buffers::copyDataRect(const hardware::Device* device, const hardware::buffers::Buffer* dest, const hardware::buffers::Buffer* orig, const size_t *dest_origin, const size_t *src_origin, const size_t *region, size_t dest_row_pitch, size_t dest_slice_pitch, size_t src_row_pitch, size_t src_slice_pitch, const std::vector<hardware::SynchronizationEvent> & events)
{
	auto const raw_events = get_raw_events(events);
	cl_event const * const raw_events_p = raw_events.size() > 0 ? raw_events.data() : nullptr;

	cl_event event_cl;
	cl_int err = clEnqueueCopyBufferRect(*device, *orig->get_cl_buffer(), *dest->get_cl_buffer(), src_origin, dest_origin, region, src_row_pitch, src_slice_pitch, dest_row_pitch, dest_slice_pitch, raw_events.size(), raw_events_p, &event_cl);
	if(err) {
		throw hardware::OpenclException(err, "clEnqueueCopyBufferRect", __FILE__, __LINE__);
	}

	const hardware::SynchronizationEvent new_event(event_cl);
	err = clReleaseEvent(event_cl);
	if(err) {
		throw hardware::OpenclException(err, "clReleaseEvent", __FILE__, __LINE__);
	}
	return new_event;
}

#ifdef CL_VERSION_1_2
void hardware::buffers::Buffer::migrate(hardware::Device * device, const std::vector<hardware::SynchronizationEvent>& events, cl_mem_migration_flags flags)
{
	size_t num_events = events.size();
	std::vector<cl_event> cl_events(num_events);
	for(size_t i = 0; i < num_events; ++i) {
		cl_events[i] = events[i].raw();
	}
	cl_event * events_p = (num_events > 0) ? &cl_events[0] : 0;

	cl_int err = clEnqueueMigrateMemObjects(*device, 1, this->get_cl_buffer(), flags, num_events, events_p, 0);
	if(err) {
		throw hardware::OpenclException(err, "clEnqueueMigrateMemoryObjects", __FILE__, __LINE__);
	}

	// update device used by this buffer
	this->device = device;
}
#endif

void memObjectReleased(cl_mem, void * user_data)
{
	MemObjectAllocationTracer * release_info = static_cast<MemObjectAllocationTracer *>(user_data);
	delete release_info;
}

std::unique_ptr<hardware::buffers::MappedBufferHandle> hardware::buffers::Buffer::map(cl_map_flags flags) const
{
	using hardware::buffers::MappedBufferHandle;

	cl_event raw_event;
	cl_int err;
	void * mapped_mem = clEnqueueMapBuffer(*device, cl_buffer, CL_FALSE, flags, 0, bytes, 0, nullptr, &raw_event, &err);
	if(err) {
		throw hardware::OpenclException(err, "clEnqueueMapBuffer", __FILE__, __LINE__);
	}
	hardware::SynchronizationEvent event(raw_event);
	err = clReleaseEvent(raw_event);
	if(err) {
		throw hardware::OpenclException(err, "clReleaseEvent", __FILE__, __LINE__);
	}
	return std::unique_ptr<MappedBufferHandle>(new MappedBufferHandle(cl_buffer, *device, mapped_mem, event));
}

hardware::buffers::MappedBufferHandle::MappedBufferHandle(cl_mem buf, cl_command_queue queue, void * mapped_ptr, hardware::SynchronizationEvent map_event)
 : buf(buf), queue(queue), mapped_ptr(mapped_ptr), map_event(map_event) { }
hardware::buffers::MappedBufferHandle::~MappedBufferHandle()
{
	cl_int err = clEnqueueUnmapMemObject(queue, buf, mapped_ptr, 0, nullptr, nullptr);
	if(err) {
		throw hardware::OpenclException(err, "clEnqueueUnmapMemObject", __FILE__, __LINE__);
	}
}
void * hardware::buffers::MappedBufferHandle::get_mapped_ptr() const
{
	return mapped_ptr;
}
hardware::SynchronizationEvent hardware::buffers::MappedBufferHandle::get_map_event() const
{
	return map_event;
}
