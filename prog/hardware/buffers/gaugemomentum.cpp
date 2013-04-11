/** @file
 * Implementation of the hardware::buffers::Gaugemomentum class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "gaugemomentum.hpp"
#include "../device.hpp"
#include "../system.hpp"
#include "../code/gaugemomentum.hpp"

#include <stdexcept>

typedef hmc_float soa_storage_t;
const size_t soa_storage_lanes = 8;

static size_t calculate_gaugemomentum_buffer_size(size_t elems, hardware::Device * device);

hardware::buffers::Gaugemomentum::Gaugemomentum(size_t elems, hardware::Device * device)
	: Buffer(calculate_gaugemomentum_buffer_size(elems, device), device),
	  elems(elems),
	  soa(check_Gaugemomentum_for_SOA(device))
{
	// nothing to do
}

size_t hardware::buffers::check_Gaugemomentum_for_SOA(hardware::Device * device)
{
	return device->get_prefers_soa();
}

static size_t calculate_gaugemomentum_buffer_size(size_t elems, hardware::Device * device)
{
	using namespace hardware::buffers;
	if(check_Gaugemomentum_for_SOA(device)) {
		size_t stride = get_Gaugemomentum_buffer_stride(elems, device);
		return stride * soa_storage_lanes * sizeof(soa_storage_t);
	} else {
		return elems * sizeof(ae);
	}
}

size_t hardware::buffers::get_Gaugemomentum_buffer_stride(size_t elems, Device * device)
{
	return device->recommend_stride(elems, sizeof(soa_storage_t), soa_storage_lanes);
}

size_t hardware::buffers::Gaugemomentum::get_elements() const noexcept
{
	return elems;
}

bool hardware::buffers::Gaugemomentum::is_soa() const noexcept
{
	return soa;
}

void hardware::buffers::Gaugemomentum::load(const ae * ptr, size_t elems, size_t offset) const
{
	if(is_soa()) {
		throw std::logic_error("Data cannot be loaded into SOA buffers.");
	} else {
		Buffer::load(ptr, elems * sizeof(ae), offset * sizeof(ae));
	}
}

void hardware::buffers::Gaugemomentum::dump(ae * ptr, size_t elems, size_t offset) const
{
	if(is_soa()) {
		auto device = get_device();
		auto gm_code = device->get_gaugemomentum_code();
		gm_code->exportGaugemomentumBuffer(ptr, this);
	} else {
		Buffer::dump(ptr, elems * sizeof(ae), offset * sizeof(ae));
	}
}

void hardware::buffers::Gaugemomentum::load_raw(const void * ptr, size_t bytes, size_t offset) const
{
	logger.trace() << "Loading raw data into Gaugemomentum buffer.";
	Buffer::load(ptr, bytes, offset);
}

void hardware::buffers::Gaugemomentum::dump_raw(void * ptr, size_t bytes, size_t offset) const
{
	logger.trace() << "Dumping raw data from Gaugemomentum buffer.";
	Buffer::dump(ptr, bytes, offset);
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
