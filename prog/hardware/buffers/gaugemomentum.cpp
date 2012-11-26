/** @file
 * Implementation of the hardware::buffers::Gaugemomentum class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "gaugemomentum.hpp"
#include "../device.hpp"
#include "../system.hpp"

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

void hardware::buffers::Gaugemomentum::load(const ae * ptr) const
{
	if(is_soa()) {
		throw std::logic_error("Data cannot be loaded into SOA buffers.");
	} else {
		Buffer::load(ptr);
	}
}

void hardware::buffers::Gaugemomentum::dump(ae * ptr) const
{
	if(is_soa()) {
		throw std::logic_error("Data cannot be dumped from SOA buffers.");
	} else {
		Buffer::dump(ptr);
	}
}
