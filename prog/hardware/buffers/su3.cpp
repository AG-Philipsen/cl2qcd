/** @file
 * Implementation of the hardware::buffers::SU3 class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "su3.hpp"
#include "../device.hpp"

#include <stdexcept>

typedef hmc_complex soa_storage_t;
const size_t soa_storage_lanes = 9;

static size_t calculate_su3_buffer_size(size_t elems, hardware::Device * device);

hardware::buffers::SU3::SU3(size_t elems, hardware::Device * device)
	: Buffer(calculate_su3_buffer_size(elems, device), device),
	  elems(elems),
	  soa(check_SU3_for_SOA(device))
{
	// nothing to do
}

size_t hardware::buffers::check_SU3_for_SOA(hardware::Device * device)
{
	return device->get_prefers_soa();
}

static size_t calculate_su3_buffer_size(size_t elems, hardware::Device * device)
{
	using namespace hardware::buffers;
	if(check_SU3_for_SOA(device)) {
		size_t stride = get_SU3_buffer_stride(elems, device);
		return stride * soa_storage_lanes * sizeof(soa_storage_t);
	} else {
		return elems * sizeof(Matrixsu3);
	}
}

size_t hardware::buffers::get_SU3_buffer_stride(size_t elems, Device * device)
{
	return device->recommend_stride(elems, sizeof(soa_storage_t), soa_storage_lanes);
}

size_t hardware::buffers::SU3::get_elements() const noexcept
{
	return elems;
}

bool hardware::buffers::SU3::is_soa() const noexcept
{
	return soa;
}

void hardware::buffers::SU3::load(const Matrixsu3 * ptr, size_t elems, size_t offset) const
{
	if(is_soa()) {
		throw std::logic_error("Data cannot be loaded into SOA buffers.");
	} else {
		Buffer::load(ptr, elems * sizeof(Matrixsu3), offset * sizeof(Matrixsu3));
	}
}

void hardware::buffers::SU3::dump(Matrixsu3 * ptr, size_t elems, size_t offset) const
{
	if(is_soa()) {
		throw std::logic_error("Data cannot be dumped from SOA buffers.");
	} else {
		Buffer::dump(ptr, elems * sizeof(Matrixsu3), offset * sizeof(Matrixsu3));
	}
}

void hardware::buffers::SU3::load_raw(const void * ptr, size_t bytes, size_t offset) const
{
	logger.trace() << "Loading raw data into SU3 buffer.";
	Buffer::load(ptr, bytes, offset);
}

void hardware::buffers::SU3::dump_raw(void * ptr, size_t bytes, size_t offset) const
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

