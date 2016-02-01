/** @file
 * Implementation of the hardware::buffers::Gaugemomentum class
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

#include "gaugemomentum.hpp"
#include "../device.hpp"
#include "../system.hpp"
#include "../code/gaugemomentum.hpp"

#include <stdexcept>

typedef hmc_float soa_storage_t;
const size_t soa_storage_lanes = 8;

static size_t calculate_gaugemomentum_buffer_size(size_t elems, const hardware::Device * device);

//redundant
size_t hardware::buffers::calculateGaugemomentumSize(LatticeExtents latticeExtentsIn) noexcept
{
	return 	latticeExtentsIn.getLatticeVolume() * NDIM;
}

hardware::buffers::Gaugemomentum::Gaugemomentum(const size_t elems, const hardware::Device * device)
	: Buffer(calculate_gaugemomentum_buffer_size(elems, device), device),
	  elems(elems),
	  soa(check_Gaugemomentum_for_SOA(device))
{}

hardware::buffers::Gaugemomentum::Gaugemomentum(const LatticeExtents lE, const hardware::Device * device)
	: Buffer(calculate_gaugemomentum_buffer_size(calculateGaugemomentumSize(lE), device), device),
	  elems(calculateGaugemomentumSize(lE)),
	  soa(check_Gaugemomentum_for_SOA(device))
{}

size_t hardware::buffers::check_Gaugemomentum_for_SOA(const hardware::Device * device)
{
	return device->get_prefers_soa();
}

static size_t calculate_gaugemomentum_buffer_size(const size_t elems, const hardware::Device * device)
{
	using namespace hardware::buffers;
	if(check_Gaugemomentum_for_SOA(device))
	{
		size_t stride = get_Gaugemomentum_buffer_stride(elems, device);
		return stride * soa_storage_lanes * sizeof(soa_storage_t);
	}
	else
	{
		return elems * sizeof(ae);
	}
}

size_t hardware::buffers::get_Gaugemomentum_buffer_stride(const size_t elems, const Device * device)
{
	return device->recommendStride(elems, sizeof(soa_storage_t), soa_storage_lanes);
}

size_t hardware::buffers::Gaugemomentum::get_elements() const noexcept
{
	return elems;
}

bool hardware::buffers::Gaugemomentum::is_soa() const noexcept
{
	return soa;
}

void hardware::buffers::Gaugemomentum::load(const ae * ptr, const size_t elems, const size_t offset) const
{
	if(is_soa()) {
		throw std::logic_error("Data cannot be loaded into SOA buffers.");
	} else {
		Buffer::load(ptr, elems * sizeof(ae), offset * sizeof(ae));
	}
}

void hardware::buffers::Gaugemomentum::dump(ae * ptr, const size_t elems, const size_t offset) const
{
	if(is_soa()) {
		auto device = get_device();
		auto gm_code = device->getGaugemomentumCode();
		gm_code->exportGaugemomentumBuffer(ptr, this);
	} else {
		Buffer::dump(ptr, elems * sizeof(ae), offset * sizeof(ae));
	}
}

void hardware::buffers::Gaugemomentum::load_raw(const void * ptr, const size_t bytes, const size_t offset) const
{
	logger.trace() << "Loading raw data into Gaugemomentum buffer.";
	Buffer::load(ptr, bytes, offset);
}

void hardware::buffers::Gaugemomentum::dump_raw(void * ptr, const size_t bytes, const size_t offset) const
{
	logger.trace() << "Dumping raw data from Gaugemomentum buffer.";
	Buffer::dump(ptr, bytes, offset);
}

void hardware::buffers::Gaugemomentum::loadRect_raw(const void* src, const size_t *buffer_origin, const size_t *host_origin, const size_t *region, size_t buffer_row_pitch, size_t buffer_slice_pitch, size_t host_row_pitch, size_t host_slice_pitch) const
{
	logger.trace() << "Loading raw data into Gaugemomentum buffer.";
	Buffer::load_rect(src, buffer_origin, host_origin, region, buffer_row_pitch, buffer_slice_pitch, host_row_pitch, host_slice_pitch);
}

void hardware::buffers::Gaugemomentum::dumpRect_raw(void* dest, const size_t *buffer_origin, const size_t *host_origin, const size_t *region, size_t buffer_row_pitch, size_t buffer_slice_pitch, size_t host_row_pitch, size_t host_slice_pitch) const
{
	logger.trace() << "Dumping raw data from Gaugemomentum buffer.";
	Buffer::dump_rect(dest, buffer_origin, host_origin, region, buffer_row_pitch, buffer_slice_pitch, host_row_pitch, host_slice_pitch);
}

hardware::SynchronizationEvent hardware::buffers::Gaugemomentum::loadRect_rawAsync(const void* src, const size_t *buffer_origin, const size_t *host_origin, const size_t *region, size_t buffer_row_pitch, size_t buffer_slice_pitch, size_t host_row_pitch, size_t host_slice_pitch, const hardware::SynchronizationEvent& event) const
{
	logger.trace() << "Loading raw data into Gaugemomentum buffer.";
	return Buffer::load_rectAsync(src, buffer_origin, host_origin, region, buffer_row_pitch, buffer_slice_pitch, host_row_pitch, host_slice_pitch, event);
}

hardware::SynchronizationEvent hardware::buffers::Gaugemomentum::dumpRect_rawAsync(void* dest, const size_t *buffer_origin, const size_t *host_origin, const size_t *region, size_t buffer_row_pitch, size_t buffer_slice_pitch, size_t host_row_pitch, size_t host_slice_pitch, const hardware::SynchronizationEvent& event) const
{
	logger.trace() << "Dumping raw data from Gaugemomentum buffer.";
	return Buffer::dump_rectAsync(dest, buffer_origin, host_origin, region, buffer_row_pitch, buffer_slice_pitch, host_row_pitch, host_slice_pitch, event);
}


size_t hardware::buffers::Gaugemomentum::get_storage_type_size() const noexcept
{
	return soa ? sizeof(soa_storage_t) : sizeof(ae);
}

size_t hardware::buffers::Gaugemomentum::get_lane_stride() const noexcept
{
	return soa ? (get_bytes() / sizeof(soa_storage_t) / soa_storage_lanes) : 0;
}

size_t hardware::buffers::Gaugemomentum::get_lane_count() const noexcept
{
	return soa ? soa_storage_lanes : 1;
}
