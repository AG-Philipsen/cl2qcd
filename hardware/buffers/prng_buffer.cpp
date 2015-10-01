/** @file
 * Implementation of the hardware::buffers::Buffer class
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

#include "prng_buffer.hpp"

hardware::buffers::PRNGBuffer::PRNGBuffer(size_t elems, Device * device)
	: hardware::buffers::Buffer{elems * sizeof(prng_state_t), device}, elems(elems) {}

hardware::buffers::PRNGBuffer::PRNGBuffer(Device * device, const bool useSameRandomNumbers)
	: hardware::buffers::PRNGBuffer::PRNGBuffer{get_prng_buffer_size(device, useSameRandomNumbers), device} {}

size_t hardware::buffers::PRNGBuffer::get_elements() const noexcept
{
	return elems;
}

void hardware::buffers::PRNGBuffer::load(const prng_state_t * array) const
{
	Buffer::load(array);
}

void hardware::buffers::PRNGBuffer::dump(prng_state_t * array) const
{
	Buffer::dump(array);
}

size_t hardware::buffers::get_prng_buffer_size(const hardware::Device * device, const bool useSameRandomNumbers)
{
	if(useSameRandomNumbers) {
		return 1.;
	} else {
#ifdef USE_PRNG_RANLUX
		// make num of random states equal to default num of global threads
		// TODO make this somewhat more automatic (avoid code duplication)
		if(device->get_device_type() == CL_DEVICE_TYPE_GPU) {
			return 4 * device->get_preferred_local_thread_num() * device->get_num_compute_units();
		} else {
			return device->get_preferred_local_thread_num() * device->get_num_compute_units();
		}
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX
	}
}

