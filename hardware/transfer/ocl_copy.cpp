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

#include "ocl_copy.hpp"

#include <stdexcept>
#include "../../exceptions.h"
#include "../device.hpp"

namespace {
	size_t get_required_buffer_size(const size_t * region);
}

hardware::transfer::OclCopy::OclCopy(hardware::Device * const from, hardware::Device * const to)
	: Transfer(from, to), transfer_buffers(), load_event(), dump_event()
{
	// nothing to do
}

hardware::transfer::OclCopy::~OclCopy()
{
	// nothing to do
}

hardware::SynchronizationEvent hardware::transfer::OclCopy::load(const hardware::buffers::Buffer* orig, const size_t *src_origin, const size_t *region, size_t src_row_pitch, size_t src_slice_pitch, const hardware::SynchronizationEvent& event)
{
	auto const transfer_buffer = get_transfer_buffer(region);

	// the transfer may neither overlap with a dump or a load, as both use the transfer buffer
	const size_t transfer_buffer_origin[] = { 0, 0, 0 };
	auto * const device = get_src_device();
	load_event = copyDataRect(device, transfer_buffer, orig, transfer_buffer_origin, src_origin, region, 0, 0, src_row_pitch, src_slice_pitch, {load_event, dump_event, event});
	device->flush();

	return load_event;
}

hardware::SynchronizationEvent hardware::transfer::OclCopy::transfer()
{
	// nothing to do
	// if the load is complet "transfer" is also complete (transfer complete means dump can start)
	return load_event;
}

hardware::SynchronizationEvent hardware::transfer::OclCopy::dump(const hardware::buffers::Buffer* dest, const size_t *dest_origin, const size_t *region, size_t dest_row_pitch, size_t dest_slice_pitch, const hardware::SynchronizationEvent& event)
{
	auto const transfer_buffer = get_transfer_buffer(region);

	// the transfer may neither overlap with a dump or a load, as both use the transfer buffer
	const size_t transfer_buffer_origin[] = { 0, 0, 0 };
	auto * const device = get_dest_device();
	dump_event = copyDataRect(device, dest, transfer_buffer, dest_origin, transfer_buffer_origin, region, dest_row_pitch, dest_slice_pitch, 0, 0, {load_event, dump_event, event});
	device->flush();

	return dump_event;
}

hardware::buffers::Buffer * hardware::transfer::OclCopy::get_transfer_buffer(const size_t * const region)
{
	auto const required_buffer_size = get_required_buffer_size(region);
	auto & transfer_buffer_handle = transfer_buffers[required_buffer_size];
	if(! (bool) transfer_buffer_handle) {
		transfer_buffer_handle = std::unique_ptr<hardware::buffers::Buffer>(new hardware::buffers::Buffer(required_buffer_size, get_src_device(), false));
	}
	return transfer_buffer_handle.get();
}

namespace {
	size_t get_required_buffer_size(const size_t * const region)
	{
		return region[0] * region[1] * region[2];
	}
}
