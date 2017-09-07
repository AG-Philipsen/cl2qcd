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

#ifndef _HARDWARE_LATTICES_PRNG_
#define _HARDWARE_LATTICES_PRNG_

#include "../system.hpp"
#include "../buffers/prng_buffer.hpp"

namespace hardware {

namespace lattices {

class PRNG
{
public:
	PRNG(const hardware::System& system, uint32_t seed, bool useSameRandomNumbers);

	virtual ~PRNG();
	PRNG& operator=(const PRNG&) = delete;
	PRNG(const PRNG&) = delete;
	PRNG() = delete;

	const std::vector<const hardware::buffers::PRNGBuffer*> get_buffers() const noexcept;

private:
	std::vector<const hardware::buffers::PRNGBuffer*> buffers;
	const hardware::System& system;
};
}

}
#endif /* _HARDWARE_LATTICES_PRNG_ */
