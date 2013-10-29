/** @file
 * Implementation of the hardware::buffers::Spinor class
 *
 * (c) 2013 Alessandro Sciarra <sciarra@compeng.uni-frankfurt.de>
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "su3vec.hpp"
#include "../device.hpp"
#include "plain.hpp"
#include "../code/spinors_staggered.hpp"

#include <stdexcept>

typedef hmc_complex soa_storage_t;
const size_t soa_storage_lanes = 3;

static size_t calculate_su3vec_buffer_size(size_t elems, hardware::Device * device);

hardware::buffers::SU3vec::SU3vec(size_t elems, hardware::Device * device)
	: Buffer(calculate_su3vec_buffer_size(elems, device), device),
	  elems(elems),
	  soa(check_su3vec_for_SOA(device))
{
	// nothing to do
}

size_t hardware::buffers::check_su3vec_for_SOA(hardware::Device * device)
{
	return device->get_prefers_soa();
}

static size_t calculate_su3vec_buffer_size(size_t elems, hardware::Device * device)
{
	using namespace hardware::buffers;
	if(check_su3vec_for_SOA(device)) {
		size_t stride = get_su3vec_buffer_stride(elems, device);
		return stride * soa_storage_lanes * sizeof(soa_storage_t);
	} else {
		return elems * sizeof(su3vec);
	}
}

size_t hardware::buffers::get_su3vec_buffer_stride(size_t elems, Device * device)
{
	return device->recommend_stride(elems, sizeof(soa_storage_t), soa_storage_lanes);
}

size_t hardware::buffers::SU3vec::get_elements() const noexcept
{
	return elems;
}

bool hardware::buffers::SU3vec::is_soa() const noexcept
{
	return soa;
}

void hardware::buffers::SU3vec::load(const su3vec * ptr, size_t elems, size_t offset) const
{
	if(is_soa()) {
		auto device = get_device();
		Plain<su3vec> plain(get_elements(), device);
		plain.load(ptr, elems * sizeof(su3vec), offset * sizeof(su3vec));
		device->get_spinor_staggered_code()->convert_staggered_field_to_SoA_eo_device(this, &plain);
		device->synchronize();
	} else {
		Buffer::load(ptr, elems * sizeof(su3vec), offset * sizeof(su3vec));
	}
}

void hardware::buffers::SU3vec::dump(su3vec * ptr, size_t elems, size_t offset) const
{
	if(is_soa()) {
		auto device = get_device();
		Plain<su3vec> plain(get_elements(), device);
		device->get_spinor_staggered_code()->convert_staggered_field_from_SoA_eo_device(&plain, this);
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
	return soa ? sizeof(soa_storage_t) : sizeof(spinor);
}

size_t hardware::buffers::SU3vec::get_lane_stride() const noexcept
{
	return soa ? (get_bytes() / sizeof(soa_storage_t) / soa_storage_lanes) : 0;
}

size_t hardware::buffers::SU3vec::get_lane_count() const noexcept
{
	return soa ? soa_storage_lanes : 1;
}
