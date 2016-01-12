/** @file
 * Implementation of DirectGMA based transfer method
 *
 * Copyright (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
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

#include "dgma.hpp"

#ifdef CL_MEM_BUS_ADDRESSABLE_AMD // make sure definitions for DGMA are available

#include <stdexcept>
#include "../../executables/exceptions.h"
#include "../device.hpp"
#include "../system.hpp"

#define DGMA_EXTENSION "cl_amd_bus_addressable_memory"

namespace {
size_t get_required_buffer_size(const size_t * region);
}

namespace hardware {
namespace transfer {

class DGMAGhostBuffer {
public:
	DGMAGhostBuffer(cl_command_queue device, hardware::buffers::Buffer * src, hardware::System const & system);
	DGMAGhostBuffer(DGMAGhostBuffer const &) = delete;
	~DGMAGhostBuffer();
	cl_mem get_cl_buffer() const;

private:
	size_t size;
	cl_mem buffer;
};

}
}

hardware::transfer::DirectGMA::DirectGMA(hardware::Device * const from, hardware::Device * const to, hardware::System const & system)
	: Transfer(from, to), src_cache(), ghost(), dest_cache(), load_event(), transfer_event(), dump_event(), active_size(0), system(system)
{
	if(!from->check_extension(DGMA_EXTENSION) || !to->check_extension(DGMA_EXTENSION)) {
		logger.error() << "DirectGMA is not supported by the used devices.";
		throw DGMAUnsupported();
	}
	cl_int err;
	transfer_queue = clCreateCommandQueue(system.getContext(), from->get_id(), 0, &err);
	if(err) {
		logger.error() << "Failed to create command queue for asynchroneous transfers. OpenCL Error: " << err;
		throw Opencl_Error(err, "clEnqueueCopyBuffer", __FILE__, __LINE__);
	}
}

hardware::transfer::DirectGMA::~DirectGMA()
{
	clReleaseCommandQueue(transfer_queue);
}

hardware::SynchronizationEvent hardware::transfer::DirectGMA::load(const hardware::buffers::Buffer* orig, const size_t *src_origin, const size_t *region, size_t src_row_pitch, size_t src_slice_pitch, const hardware::SynchronizationEvent& event)
{
	active_size = get_required_buffer_size(region);
	ensure_buffers(active_size);

	// the transfer may neither overlap with a dump or a load, as both use the transfer buffer
	const size_t transfer_buffer_origin[] = { 0, 0, 0 };
	auto * const device = get_src_device();
	load_event = copyDataRect(device, src_cache.get(), orig, transfer_buffer_origin, src_origin, region, 0, 0, src_row_pitch, src_slice_pitch, {load_event, transfer_event, event});
	device->flush(); // event and buffer will be used by other queues

	return load_event;
}

hardware::SynchronizationEvent hardware::transfer::DirectGMA::transfer()
{
	if(active_size == 0) {
		logger.error() << "Transfer requested without first loading data into the link.";
		throw std::logic_error("Transfer requested without first loading data into the link.");
	}

	// transfer data to destination device
	auto const events = get_raw_events( {load_event, transfer_event, dump_event});
	cl_event const * const events_p = (events.size() > 0) ? events.data() : nullptr;

	cl_event raw_transfer_event;
	cl_int clerr = clEnqueueCopyBuffer(transfer_queue, *src_cache->get_cl_buffer(), ghost->get_cl_buffer(), 0, 0, active_size, events.size(), events_p, &raw_transfer_event);
	if(clerr) {
		logger.error() << "Failed to transfer data between devices using seperate transfer queue. OpenCL Error: " << clerr;
		throw Opencl_Error(clerr, "clEnqueueCopyBuffer", __FILE__, __LINE__);
	}

	transfer_event = SynchronizationEvent(raw_transfer_event);
	clerr = clReleaseEvent(raw_transfer_event);
	if(clerr) {
		throw Opencl_Error(clerr, "clReleaseEvent", __FILE__, __LINE__);
	}

	// event and buffer will be used by other queues -> flush
	clerr = clFlush(transfer_queue);
	if(clerr) {
		logger.error() << "Failed to transfer data between devices using seperate transfer queue. OpenCL Error: " << clerr;
		throw Opencl_Error(clerr, "clEnqueueCopyBuffer", __FILE__, __LINE__);
	}

	return transfer_event;
}

hardware::SynchronizationEvent hardware::transfer::DirectGMA::dump(const hardware::buffers::Buffer* dest, const size_t *dest_origin, const size_t *region, size_t dest_row_pitch, size_t dest_slice_pitch, const hardware::SynchronizationEvent& event)
{
	auto const required_buffer_size = get_required_buffer_size(region);
	if(required_buffer_size != active_size) {
		logger.error() << "Buffer requested to be dumped has a different size than the last buffer loaded";
		throw std::logic_error("Buffer requested to be dumped has a different size than the last buffer loaded");
	}

	// the transfer may neither overlap with a dump or a load, as both use the transfer buffer
	const size_t transfer_buffer_origin[] = { 0, 0, 0 };
	auto * const device = get_dest_device();
	dump_event = copyDataRect(device, dest, dest_cache.get(), dest_origin, transfer_buffer_origin, region, dest_row_pitch, dest_slice_pitch, 0, 0, {transfer_event, dump_event, event});
	device->flush(); // event might be used by other queues

	return dump_event;
}

void hardware::transfer::DirectGMA::ensure_buffers(size_t buffer_size)
{
	if(! (bool) src_cache || src_cache->get_bytes() < buffer_size) {
		src_cache = std::unique_ptr<hardware::buffers::Buffer>(new hardware::buffers::Buffer(buffer_size, get_src_device(), false));
		dest_cache = std::unique_ptr<hardware::buffers::Buffer>(new hardware::buffers::Buffer(buffer_size, get_dest_device(), false, CL_MEM_BUS_ADDRESSABLE_AMD));
		ghost = std::unique_ptr<hardware::transfer::DGMAGhostBuffer>(new hardware::transfer::DGMAGhostBuffer(transfer_queue, dest_cache.get(), system));
	}
}

hardware::transfer::DGMAGhostBuffer::DGMAGhostBuffer(cl_command_queue device, hardware::buffers::Buffer * src, hardware::System const & system)
{
	size = src->get_bytes();

	cl_bus_address_amd address;


	clEnqueueMakeBuffersResidentAMD_fn clEnqueueMakeBuffersResidentAMD = reinterpret_cast<clEnqueueMakeBuffersResidentAMD_fn>(clGetExtensionFunctionAddressForPlatform(system.get_platform(), "clEnqueueMakeBuffersResidentAMD"));
	if(!clEnqueueMakeBuffersResidentAMD) {
		throw std::runtime_error("Failed to resolve DirectGMA functions. Does the current platform support DirectGMA?");
	}

	// make sure source buffer is physically addressable
	cl_mem tmp = *src->get_cl_buffer();
	cl_int err = clEnqueueMakeBuffersResidentAMD(src->get_device()->command_queue, 1, &tmp, CL_TRUE, &address, 0, nullptr, nullptr);
	if(err) {
		throw Opencl_Error(err, "clEnqueueMakeBuffersResidentAMD", __FILE__, __LINE__);
	}

	// create the shadow buffer on the source device that is backed by the bus_addressable_buffer's memory
	buffer = clCreateBuffer(system.getContext(), CL_MEM_EXTERNAL_PHYSICAL_AMD, size, &address, &err);
	if(err) {
		logger.error() << "Failed to create ghost buffer";
		throw Opencl_Error(err, "clCreateBuffer", __FILE__, __LINE__);
	}
	err = clEnqueueMigrateMemObjects(device, 1, &buffer, 0, 0, nullptr, nullptr);
	if(err) {
		throw Opencl_Error(err, "clEnqueueMigrateMemObjects", __FILE__, __LINE__);
	}
	err = clFinish(device);
	if(err) {
		throw Opencl_Error(err, "clFinish", __FILE__, __LINE__);
	}
}

hardware::transfer::DGMAGhostBuffer::~DGMAGhostBuffer()
{
	clReleaseMemObject(buffer);
}

cl_mem hardware::transfer::DGMAGhostBuffer::get_cl_buffer() const
{
	return buffer;
}

namespace {
size_t get_required_buffer_size(const size_t * const region)
{
	return region[0] * region[1] * region[2];
}
}

#endif
