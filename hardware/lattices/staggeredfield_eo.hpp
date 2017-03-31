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

#ifndef _HARDWARE_LATTICES_STAGGEREDFIELD_EO_
#define _HARDWARE_LATTICES_STAGGEREDFIELD_EO_

#include "../system.hpp"
#include "../buffers/su3vec.hpp"
#include "../../common_header_files/types_fermions.h"

namespace hardware {

namespace lattices {

class Staggeredfield_eo {

public:
	Staggeredfield_eo(const hardware::System&);

	virtual ~Staggeredfield_eo();

    /**
     * Staggeredfield_eo cannot be copied but only move-constructed
     */
    Staggeredfield_eo() = delete;
    Staggeredfield_eo(const Staggeredfield_eo&) = delete;
    Staggeredfield_eo& operator=(const Staggeredfield_eo&) = delete;
    Staggeredfield_eo(Staggeredfield_eo&&); //own implementation since class own resources
    Staggeredfield_eo& operator=(Staggeredfield_eo&&) = delete;

	std::vector<const hardware::buffers::SU3vec *> allocate_buffers();

	const std::vector<const hardware::buffers::SU3vec *> get_buffers() const noexcept;

	void update_halo() const;

	void import(const su3vec * const host) const;

private:
	hardware::System const& system;
	const std::vector<const hardware::buffers::SU3vec *> buffers;
};
	
}
}
#endif /*_HARDWARE_LATTICES_SPINORFIELD_EO_ */
