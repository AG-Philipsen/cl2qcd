/** @file
 * Implementation of the simple opencl transfer method
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

#include "async_ocl_copy.hpp"

#include <stdexcept>
#include "../../executables/exceptions.h"
#include "../device.hpp"
#include "../system.hpp"

namespace {
size_t get_required_buffer_size(const size_t * region);
}

hardware::transfer::AsyncOclCopy::AsyncOclCopy(hardware::Device * const from, hardware::Device * const to, hardware::System const & system)
	: Transfer(from, to), src_cache(), dest_cache(), load_event(), transfer_event(), dump_event(), back_migration_event(), active_size(0)
{
	cl_int err;
	transfer_queue = clCreateCommandQueue(system.getContext(), to->get_id(), 0, &err);
	if(err) {
		logger.error() << "Failed to create command queue for asynchroneous transfers. OpenCL Error: " << err;
		throw Opencl_Error(err, "clEnqueueCopyBuffer", __FILE__, __LINE__);
	}
	back_migration_queue = clCreateCommandQueue(system.getContext(), from->get_id(), 0, &err);
	if(err) {
		logger.error() << "Failed to create command queue for asynchroneous transfers. OpenCL Error: " << err;
		throw Opencl_Error(err, "clEnqueueCopyBuffer", __FILE__, __LINE__);
	}
}

hardware::transfer::AsyncOclCopy::~AsyncOclCopy()
{
	clReleaseCommandQueue(transfer_queue);
	clReleaseCommandQueue(back_migration_queue);
}

hardware::SynchronizationEvent hardware::transfer::AsyncOclCopy::load(const hardware::buffers::Buffer* orig, const size_t *src_origin, const size_t *region, size_t src_row_pitch, size_t src_slice_pitch, const hardware::SynchronizationEvent& event)
{
	active_size = get_required_buffer_size(region);
	auto const transfer_buffer = get_src_cache(active_size);

	// the transfer may neither overlap with a dump or a load, as both use the transfer buffer
	const size_t transfer_buffer_origin[] = { 0, 0, 0 };
	load_event = copyDataRect(get_src_device(), transfer_buffer, orig, transfer_buffer_origin, src_origin, region, 0, 0, src_row_pitch, src_slice_pitch, {load_event, transfer_event, event, back_migration_event});

	return load_event;
}

hardware::SynchronizationEvent hardware::transfer::AsyncOclCopy::transfer()
{
	auto const src_cache = get_src_cache(active_size);
	auto const dest_cache = get_dest_cache(active_size);

	// transfer data to destination device
	{
		auto const events = get_raw_events( {load_event, transfer_event, dump_event, back_migration_event});
		cl_event const * const events_p = (events.size() > 0) ? events.data() : nullptr;

		cl_event raw_transfer_event;
		cl_int clerr = clEnqueueCopyBuffer(transfer_queue, *src_cache->get_cl_buffer(), *dest_cache->get_cl_buffer(), 0, 0, active_size, events.size(), events_p, &raw_transfer_event);
		if(clerr) {
			logger.error() << "Failed to transfer data between devices using seperate transfer queue. OpenCL Error: " << clerr;
			throw Opencl_Error(clerr, "clEnqueueCopyBuffer", __FILE__, __LINE__);
		}

		transfer_event = SynchronizationEvent(raw_transfer_event);
		clerr = clReleaseEvent(raw_transfer_event);
		if(clerr) {
			throw Opencl_Error(clerr, "clReleaseEvent", __FILE__, __LINE__);
		}
	}

	// migrate source cache back to source device
#if CL_VERSION_1_2
	{
		auto const events = get_raw_events( {load_event, transfer_event, back_migration_event});
		cl_event const * const events_p = (events.size() > 0) ? events.data() : nullptr;

		cl_event raw_back_migration_event;
		cl_int clerr = clEnqueueMigrateMemObjects(back_migration_queue, 1, src_cache->get_cl_buffer(), CL_MIGRATE_MEM_OBJECT_CONTENT_UNDEFINED, events.size(), events_p, &raw_back_migration_event);
		if(clerr) {
			throw Opencl_Error(clerr, "clEnqueueMigrateMemObjects", __FILE__, __LINE__);
		}
		back_migration_event = SynchronizationEvent(raw_back_migration_event);
		clerr = clReleaseEvent(raw_back_migration_event);
		if(clerr) {
			throw Opencl_Error(clerr, "clReleaseEvent", __FILE__, __LINE__);
		}
	}
#else
	// back migration is a performance optimization that is not functionally required.
	// if it is not available just work
	back_migration_event = transfer_event;
#pragma message "clEnqueueMigrateMemObjects is not available on your system. This might negativly affect multi-device performance. Please upgrade to OpenCL 1.2 or higher."
#endif

	return transfer_event;
}

hardware::SynchronizationEvent hardware::transfer::AsyncOclCopy::dump(const hardware::buffers::Buffer* dest, const size_t *dest_origin, const size_t *region, size_t dest_row_pitch, size_t dest_slice_pitch, const hardware::SynchronizationEvent& event)
{
	auto const required_buffer_size = get_required_buffer_size(region);
	if(required_buffer_size != active_size) {
		logger.error() << "Buffer requested to be dumped has a different size than the last buffer loaded";
		throw std::logic_error("Buffer requested to be dumped has a different size than the last buffer loaded");
	}
	auto const transfer_buffer = get_dest_cache(required_buffer_size);

	// the transfer may neither overlap with a dump or a load, as both use the transfer buffer
	const size_t transfer_buffer_origin[] = { 0, 0, 0 };
	dump_event = copyDataRect(get_dest_device(), dest, transfer_buffer, dest_origin, transfer_buffer_origin, region, dest_row_pitch, dest_slice_pitch, 0, 0, {transfer_event, dump_event, event});

	return dump_event;
}

hardware::buffers::Buffer * hardware::transfer::AsyncOclCopy::get_src_cache(size_t buffer_size)
{
	auto & transfer_buffer_handle = src_cache[buffer_size];
	if(! (bool) transfer_buffer_handle) {
		transfer_buffer_handle = std::unique_ptr<hardware::buffers::Buffer>(new hardware::buffers::Buffer(buffer_size, get_src_device(), false));
	}
	return transfer_buffer_handle.get();
}

hardware::buffers::Buffer * hardware::transfer::AsyncOclCopy::get_dest_cache(size_t const buffer_size)
{
	auto & transfer_buffer_handle = dest_cache[buffer_size];
	if(! (bool) transfer_buffer_handle) {
		transfer_buffer_handle = std::unique_ptr<hardware::buffers::Buffer>(new hardware::buffers::Buffer(buffer_size, get_dest_device(), false));
	}
	return transfer_buffer_handle.get();
}

namespace {
size_t get_required_buffer_size(const size_t * const region)
{
	return region[0] * region[1] * region[2];
}
}
