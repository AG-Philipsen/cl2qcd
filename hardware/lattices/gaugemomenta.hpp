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

#ifndef _HARDWARE_LATTICES_GAUGEMOMENTA_
#define _HARDWARE_LATTICES_GAUGEMOMENTA_

#include "../system.hpp"
#include "../buffers/gaugemomentum.hpp"

namespace hardware {

namespace lattices {

class Gaugemomenta
{
public:

	virtual ~Gaugemomenta();
	
	Gaugemomenta(const hardware::System& system);

	Gaugemomenta& operator=(const Gaugemomenta&) = delete;
	Gaugemomenta(const Gaugemomenta&) = delete;
	Gaugemomenta() = delete;

	std::vector<const hardware::buffers::Gaugemomentum *> allocate_buffers() const;
	const std::vector<const hardware::buffers::Gaugemomentum *> get_buffers() const noexcept;

	void import(const ae * const host) const;
	void update_halo() const;
	void update_halo_soa(const std::vector<const hardware::buffers::Gaugemomentum *> buffers, const hardware::System& system) const;
	void update_halo_aos(const std::vector<const hardware::buffers::Gaugemomentum *> buffers, const hardware::System& system) const;
private:
	hardware::System const& system;
	const std::vector<const hardware::buffers::Gaugemomentum *> buffers;
};

}

}

#endif /* _HARDWARE_LATTICES_GAUGEMOMENTA_ */
