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

#ifndef _HARDWARE_LATTICES_SPINORFIELD_
#define _HARDWARE_LATTICES_SPINORFIELD_

#include "../system.hpp"
#include "../device.hpp"
#include "../buffers/plain.hpp"
#include "../../common_header_files/types_fermions.h"


namespace hardware {

namespace lattices {

class Spinorfield
{
public:

	Spinorfield(const hardware::System&, const bool place_on_host = false);

    Spinorfield& operator=(const Spinorfield&) = delete;
    Spinorfield(const Spinorfield&) = delete;
    Spinorfield() = delete;

	virtual ~Spinorfield();

	const std::vector<const hardware::buffers::Plain<spinor> *> get_buffers() const noexcept;
	std::vector<const hardware::buffers::Plain<spinor> *> allocate_buffers();
	void import(const spinor * const host) const;
	void update_halo() const;

// protected:

	void fill_buffers();
	void clear_buffers();

private:
        hardware::System const& system;
        std::vector<const hardware::buffers::Plain<spinor> *> buffers;
        const bool place_on_host;
};

}

}

#endif /* _HARDWARE_LATTICES_SPINORFIELD_ */
