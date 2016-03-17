/*
 * Copyright 2016 Francesca Cuteri
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

#include "prng.hpp"
#include "../../host_functionality/host_random.h"
#include "../code/prng.hpp"

hardware::lattices::PRNG::PRNG(const hardware::System& system, uint32_t seed, bool useSameRandomNumbers) :
	system(system)
{
	using hardware::buffers::PRNGBuffer;

	// initialize host prng
	prng_init(seed);

	// initialize devices
	for(hardware::Device * device : system.get_devices()) {
		// create a buffer for each device
		const PRNGBuffer * buffer = new PRNGBuffer(device, useSameRandomNumbers);
		auto code = device->getPrngCode();
		code->initialize(buffer, ++seed);
		buffers.push_back(buffer);
	}
}

hardware::lattices::PRNG::~PRNG()
{
	for(const hardware::buffers::PRNGBuffer * buffer : buffers)
	{
		delete buffer;
	}
}

const std::vector<const hardware::buffers::PRNGBuffer*> hardware::lattices::PRNG::get_buffers() const noexcept
{
	return buffers;
}
