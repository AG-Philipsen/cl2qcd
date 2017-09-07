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

#ifndef _HARDWARE_LATTICES_GAUGEFIELD_
#define _HARDWARE_LATTICES_GAUGEFIELD_

#include "../system.hpp"
#include "../buffers/su3.hpp"
#include "../buffers/plain.hpp"
#include "prng.hpp"


namespace hardware {

namespace lattices {

class Gaugefield 
{
public:

	virtual ~Gaugefield();

	Gaugefield(const hardware::System& system);

	Gaugefield& operator=(const Gaugefield&) = delete;
	Gaugefield(const Gaugefield&) = delete;
	Gaugefield() = delete;

	const std::vector<const hardware::buffers::SU3 *> get_buffers() const noexcept;
	std::vector<const hardware::buffers::SU3 *> allocate_buffers();
	void release_buffers(std::vector<const hardware::buffers::SU3 *>* buffers);
	void send_gaugefield_to_buffers(const Matrixsu3 * const gf_host);
	void fetch_gaugefield_from_buffers( Matrixsu3 * const gf_host);

	void update_halo() const;
	
	void set_cold() const;
	void set_hot() const;

	void smear(unsigned int smearingSteps);
	void unsmear();

private:
	hardware::System const& system;
	std::vector<const hardware::buffers::SU3 *> buffers;
	std::vector<const hardware::buffers::SU3 *> unsmeared_buffers;

	void update_halo_soa(std::vector<const hardware::buffers::SU3 *> buffers, const hardware::System& system) const;
	void update_halo_aos(std::vector<const hardware::buffers::SU3 *> buffers, const hardware::System& system) const;

	void set_cold(Matrixsu3 * field, size_t elems) const;
	void set_hot(Matrixsu3 * field, size_t elems) const;
};

}

}

#endif /* _HARDWARE_LATTICES_GAUGEFIELD_ */
