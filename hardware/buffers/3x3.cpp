/** @file
 * Implementation of the hardware::buffers::Matrix3x3 class
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

#include "3x3.hpp"
#include "../device.hpp"

#include <stdexcept>

typedef hmc_complex soa_storage_t;
const size_t soa_storage_lanes = 9;

static size_t calculate_3x3_buffer_size(const size_t elems, const hardware::Device * device);

hardware::buffers::Matrix3x3::Matrix3x3(const size_t elems, const hardware::Device * device)
	: Buffer(calculate_3x3_buffer_size(elems, device), device),
	  elems(elems),
	  soa(check_Matrix3x3_for_SOA(device))
{
	// nothing to do
}

size_t hardware::buffers::check_Matrix3x3_for_SOA(const hardware::Device * device)
{
	return device->get_prefers_soa();
}

static size_t calculate_3x3_buffer_size(const size_t elems, const hardware::Device * device)
{
	using namespace hardware::buffers;
	if(check_Matrix3x3_for_SOA(device)) {
		size_t stride = get_Matrix3x3_buffer_stride(elems, device);
		return stride * soa_storage_lanes * sizeof(soa_storage_t);
	} else {
		return elems * sizeof(::Matrix3x3);
	}
}

size_t hardware::buffers::get_Matrix3x3_buffer_stride(const size_t elems, const Device * device)
{
	return device->recommendStride(elems, sizeof(soa_storage_t), soa_storage_lanes);
}

size_t hardware::buffers::Matrix3x3::get_elements() const noexcept
{
	return elems;
}

bool hardware::buffers::Matrix3x3::is_soa() const noexcept
{
	return soa;
}
