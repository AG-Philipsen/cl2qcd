/** @file
 * Implementation of the hardware::buffers::SU3 class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "su3.hpp"

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

void hardware::buffers::SU3::load(const Matrixsu3 * ptr) const
{
	if(is_soa()) {
		throw std::logic_error("Data cannot be loaded into SOA buffers.");
	} else {
		Buffer::load(ptr);
	}
}

void hardware::buffers::SU3::dump(Matrixsu3 * ptr) const
{
	if(is_soa()) {
		throw std::logic_error("Data cannot be dumped from SOA buffers.");
	} else {
		Buffer::dump(ptr);
	}
}
