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

#include "../../hardware/system.hpp"
#include "../buffers/su3.hpp"
#include "../buffers/plain.hpp"


namespace hardware {

namespace lattices {

class Gaugefield 
{
public:
	friend hardware::Device;

	virtual ~Gaugefield();
	
	Gaugefield(const hardware::System& system);

	std::vector<const hardware::buffers::SU3 *> allocate_buffers() const;
	void release_buffers(std::vector<const hardware::buffers::SU3 *>* buffers) const;
	void send_gaugefield_to_buffers(const std::vector<const hardware::buffers::SU3 *> buffers, const Matrixsu3 * const gf_host) const;
	void fetch_gaugefield_from_buffers(const std::vector<const hardware::buffers::SU3 *> buffers, Matrixsu3 * const gf_host) const;
	void update_halo_soa(std::vector<const hardware::buffers::SU3 *> buffers, const hardware::System& system) const;
	void update_halo_aos(std::vector<const hardware::buffers::SU3 *> buffers, const hardware::System& system) const;
private:
	hardware::System const& system;
	std::vector<const hardware::buffers::SU3 *> buffers;
	std::vector<const hardware::buffers::SU3 *> unsmeared_buffers;
};

}

}

#endif /* _HARDWARE_LATTICES_GAUGEFIELD_ */