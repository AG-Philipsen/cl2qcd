/** @file
 * Implementation of the hardware::buffers::SU3 class
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

#include "su3.hpp"
#include "../device.hpp"
#include "../code/gaugefield.hpp"

#include <stdexcept>

typedef hmc_complex soa_storage_t;
const size_t soa_storage_lanes = 9;

static size_t calculate_su3_buffer_size(size_t elems, const hardware::Device * device);

size_t hardware::buffers::calculateGaugefieldSize(const LatticeExtents latticeExtentsIn) noexcept
{
	return 	latticeExtentsIn.getLatticeVolume() * NDIM;
}

hardware::buffers::SU3::SU3(const size_t elems, const hardware::Device * device)
	: Buffer(calculate_su3_buffer_size(elems, device), device),
	  elems(elems),
	  soa(check_SU3_for_SOA(device))
{
	// nothing to do
}

hardware::buffers::SU3::SU3(const LatticeExtents lE, const hardware::Device * device)
	: Buffer(calculate_su3_buffer_size(calculateGaugefieldSize(lE), device), device),
	  elems(calculateGaugefieldSize(lE)),
	  soa(check_SU3_for_SOA(device))
{
	// nothing to do
}

size_t hardware::buffers::check_SU3_for_SOA(const hardware::Device * device)
{
	return device->get_prefers_soa();
}

static size_t calculate_su3_buffer_size(const size_t elems, const hardware::Device * device)
{
	using namespace hardware::buffers;
	if(check_SU3_for_SOA(device)) {
		size_t stride = get_SU3_buffer_stride(elems, device);
		return stride * soa_storage_lanes * sizeof(soa_storage_t);
	} else {
		return elems * sizeof(Matrixsu3);
	}
}

size_t hardware::buffers::get_SU3_buffer_stride(const size_t elems, const Device * device)
{
	return device->recommendStride(elems, sizeof(soa_storage_t), soa_storage_lanes);
}

size_t hardware::buffers::SU3::get_elements() const noexcept
{
	return elems;
}

bool hardware::buffers::SU3::is_soa() const noexcept
{
	return soa;
}

void hardware::buffers::SU3::load(const Matrixsu3 * ptr, const size_t elems, const size_t offset) const
{
	if(is_soa()) {
		throw std::logic_error("Data cannot be loaded into SOA buffers.");
	} else {
		Buffer::load(ptr, elems * sizeof(Matrixsu3), offset * sizeof(Matrixsu3));
	}
}

void hardware::buffers::SU3::dump(Matrixsu3 * ptr, const size_t elems, const size_t offset) const
{
	if(is_soa()) {
		throw std::logic_error("Data cannot be dumped from SOA buffers.");
	} else {
		Buffer::dump(ptr, elems * sizeof(Matrixsu3), offset * sizeof(Matrixsu3));
	}
}

void hardware::buffers::SU3::load_raw(const void * ptr, const size_t bytes, const size_t offset) const
{
	logger.trace() << "Loading raw data into SU3 buffer.";
	Buffer::load(ptr, bytes, offset);
}

void hardware::buffers::SU3::dump_raw(void * ptr, const size_t bytes, const size_t offset) const
{
	logger.trace() << "Dumping raw data from SU3 buffer.";
	Buffer::dump(ptr, bytes, offset);
}

void hardware::buffers::SU3::loadRect_raw(const void* src, const size_t *buffer_origin, const size_t *host_origin, const size_t *region, size_t buffer_row_pitch, size_t buffer_slice_pitch, size_t host_row_pitch, size_t host_slice_pitch) const
{
	logger.trace() << "Loading raw data into SU3 buffer.";
	Buffer::load_rect(src, buffer_origin, host_origin, region, buffer_row_pitch, buffer_slice_pitch, host_row_pitch, host_slice_pitch);
}

void hardware::buffers::SU3::dumpRect_raw(void* dest, const size_t *buffer_origin, const size_t *host_origin, const size_t *region, size_t buffer_row_pitch, size_t buffer_slice_pitch, size_t host_row_pitch, size_t host_slice_pitch) const
{
	logger.trace() << "Dumping raw data from SU3 buffer.";
	Buffer::dump_rect(dest, buffer_origin, host_origin, region, buffer_row_pitch, buffer_slice_pitch, host_row_pitch, host_slice_pitch);
}
hardware::SynchronizationEvent hardware::buffers::SU3::loadRect_rawAsync(const void* src, const size_t *buffer_origin, const size_t *host_origin, const size_t *region, size_t buffer_row_pitch, size_t buffer_slice_pitch, size_t host_row_pitch, size_t host_slice_pitch, const hardware::SynchronizationEvent& event) const
{
	logger.trace() << "Loading raw data into SU3 buffer.";
	return Buffer::load_rectAsync(src, buffer_origin, host_origin, region, buffer_row_pitch, buffer_slice_pitch, host_row_pitch, host_slice_pitch, event);
}

hardware::SynchronizationEvent hardware::buffers::SU3::dumpRect_rawAsync(void* dest, const size_t *buffer_origin, const size_t *host_origin, const size_t *region, size_t buffer_row_pitch, size_t buffer_slice_pitch, size_t host_row_pitch, size_t host_slice_pitch, const hardware::SynchronizationEvent& event) const
{
	logger.trace() << "Dumping raw data from SU3 buffer.";
	return Buffer::dump_rectAsync(dest, buffer_origin, host_origin, region, buffer_row_pitch, buffer_slice_pitch, host_row_pitch, host_slice_pitch, event);
}



size_t hardware::buffers::SU3::get_storage_type_size() const noexcept
{
	return soa ? sizeof(soa_storage_t) : sizeof(Matrixsu3);
}

size_t hardware::buffers::SU3::get_lane_stride() const noexcept
{
	return soa ? (get_bytes() / sizeof(soa_storage_t) / soa_storage_lanes) : 0;
}

size_t hardware::buffers::SU3::get_lane_count() const noexcept
{
	return soa ? soa_storage_lanes : 1;
}

