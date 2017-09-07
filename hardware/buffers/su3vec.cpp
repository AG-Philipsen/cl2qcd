/** @file
 * Implementation of the hardware::buffers::Spinor class
 *
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.phys.uni-frankfurt.de>
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

#include "su3vec.hpp"
#include "../device.hpp"
#include "plain.hpp"
#include "../code/spinors_staggered.hpp"

#include <stdexcept>

typedef hmc_complex soa_storage_t;
const size_t soa_storage_lanes = 3;

static size_t calculate_su3vec_buffer_size(const size_t elems, const hardware::Device * device);

size_t hardware::buffers::calculateEvenOddSpinorfieldSize(LatticeExtents latticeExtendsIn) noexcept
{
	return latticeExtendsIn.getLatticeVolume() / 2;
}


hardware::buffers::SU3vec::SU3vec(const size_t elems, const hardware::Device * device)
	: Buffer(calculate_su3vec_buffer_size(elems, device), device),
	  elems(elems),
	  soa(check_su3vec_for_SOA(device))
{
	// nothing to do
}


hardware::buffers::SU3vec::SU3vec(const LatticeExtents lE, const hardware::Device * device, const int fact)
	: Buffer(calculate_su3vec_buffer_size(calculateEvenOddSpinorfieldSize(lE)*fact, device), device),
	  elems(calculateEvenOddSpinorfieldSize(lE)*fact),
	  soa(check_su3vec_for_SOA(device))
{
	// nothing to do
}


size_t hardware::buffers::check_su3vec_for_SOA(const hardware::Device * device)
{
	return device->get_prefers_soa();
}

static size_t calculate_su3vec_buffer_size(const size_t elems, const hardware::Device * device)
{
	using namespace hardware::buffers;
	if(check_su3vec_for_SOA(device)) {
		size_t stride = get_su3vec_buffer_stride(elems, device);
		return stride * soa_storage_lanes * sizeof(soa_storage_t);
	} else {
		return elems * sizeof(su3vec);
	}
}

size_t hardware::buffers::get_su3vec_buffer_stride(const size_t elems, const Device * device)
{
	return device->recommendStride(elems, sizeof(soa_storage_t), soa_storage_lanes);
}

size_t hardware::buffers::SU3vec::get_elements() const noexcept
{
	return elems;
}

bool hardware::buffers::SU3vec::is_soa() const noexcept
{
	return soa;
}

void hardware::buffers::SU3vec::load(const su3vec * ptr, const size_t elems, const size_t offset) const
{
	if(is_soa()) {
		auto device = get_device();
		Plain<su3vec> plain(get_elements(), device);
		plain.load(ptr, elems * sizeof(su3vec), offset * sizeof(su3vec));
		device->getSpinorStaggeredCode()->convert_staggered_field_to_SoA_eo_device(this, &plain);
		device->synchronize();
	} else {
		Buffer::load(ptr, elems * sizeof(su3vec), offset * sizeof(su3vec));
	}
}

void hardware::buffers::SU3vec::dump(su3vec * ptr, const size_t elems, const size_t offset) const
{
	if(is_soa()) {
		auto device = get_device();
		Plain<su3vec> plain(get_elements(), device);
		device->getSpinorStaggeredCode()->convert_staggered_field_from_SoA_eo_device(&plain, this);
		plain.dump(ptr, elems * sizeof(su3vec), offset * sizeof(su3vec));
	} else {
		Buffer::dump(ptr, elems * sizeof(su3vec), offset * sizeof(su3vec));
	}
}



/*To be added...

void hardware::buffers::Spinor::load_raw(const void * ptr, size_t bytes, size_t offset) const
{
	logger.trace() << "Loading raw data into Spinor buffer.";
	Buffer::load(ptr, bytes, offset);
}

void hardware::buffers::Spinor::dump_raw(void * ptr, size_t bytes, size_t offset) const
{
	logger.trace() << "Dumping raw data from Spinor buffer.";
	Buffer::dump(ptr, bytes, offset);
}

void hardware::buffers::Spinor::loadRect_raw(const void* src, const size_t *buffer_origin, const size_t *host_origin, const size_t *region, size_t buffer_row_pitch, size_t buffer_slice_pitch, size_t host_row_pitch, size_t host_slice_pitch) const
{
	logger.trace() << "Loading raw data into Spinor buffer.";
	Buffer::load_rect(src, buffer_origin, host_origin, region, buffer_row_pitch, buffer_slice_pitch, host_row_pitch, host_slice_pitch);
}

void hardware::buffers::Spinor::dumpRect_raw(void* dest, const size_t *buffer_origin, const size_t *host_origin, const size_t *region, size_t buffer_row_pitch, size_t buffer_slice_pitch, size_t host_row_pitch, size_t host_slice_pitch) const
{
	logger.trace() << "Dumping raw data from Spinor buffer.";
	Buffer::dump_rect(dest, buffer_origin, host_origin, region, buffer_row_pitch, buffer_slice_pitch, host_row_pitch, host_slice_pitch);
}

hardware::SynchronizationEvent hardware::buffers::Spinor::loadRect_rawAsync(const void* src, const size_t *buffer_origin, const size_t *host_origin, const size_t *region, size_t buffer_row_pitch, size_t buffer_slice_pitch, size_t host_row_pitch, size_t host_slice_pitch, const hardware::SynchronizationEvent& event) const
{
	logger.trace() << "Loading raw data into Spinor buffer.";
	return Buffer::load_rectAsync(src, buffer_origin, host_origin, region, buffer_row_pitch, buffer_slice_pitch, host_row_pitch, host_slice_pitch, event);
}

hardware::SynchronizationEvent hardware::buffers::Spinor::dumpRect_rawAsync(void* dest, const size_t *buffer_origin, const size_t *host_origin, const size_t *region, size_t buffer_row_pitch, size_t buffer_slice_pitch, size_t host_row_pitch, size_t host_slice_pitch, const hardware::SynchronizationEvent& event) const
{
	logger.trace() << "Dumping raw data from Spinor buffer.";
	return Buffer::dump_rectAsync(dest, buffer_origin, host_origin, region, buffer_row_pitch, buffer_slice_pitch, host_row_pitch, host_slice_pitch, event);
}
*/

size_t hardware::buffers::SU3vec::get_storage_type_size() const noexcept
{
	return soa ? sizeof(soa_storage_t) : sizeof(su3vec);
}

size_t hardware::buffers::SU3vec::get_lane_stride() const noexcept
{
	return soa ? (get_bytes() / sizeof(soa_storage_t) / soa_storage_lanes) : 0;
}

size_t hardware::buffers::SU3vec::get_lane_count() const noexcept
{
	return soa ? soa_storage_lanes : 1;
}
