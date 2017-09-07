/** @file
 * Utility functions for physics::lattices
 *
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
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

#ifndef _PHYSICS_LATTICES_UTIL_
#define _PHYSICS_LATTICES_UTIL_

#include <stdexcept>
#include "../../hardware/buffers/buffer.hpp"

namespace physics {

namespace lattices {

/**
 * Copy the contents of one lattice to another
 *
 * \param[out] dest The lattice to copy to
 * \param[in]  from The lattice to copy from
 */
template<class T> void copyData(const T* to, const T* from);

/**
 * Copy the contents of one lattice to another
 *
 * \param[out] dest The lattice to copy to
 * \param[in]  from The lattice to copy from
 */
template<class T> void copyData(const T* to, const T& from);


/**
 * Pseudo-Randomize a lattice. This can be usefull for testcases.
 *
 * @todo implement for multi-buffer
 */
template <class Lattice, typename Basetype> void pseudo_randomize(const Lattice* to, int seed);



/*
 * TEMPLATE IMPLEMENTATIONS
 */
template<class T> void copyData(const T* to, const T& from)
{
	auto from_buffers = from.get_buffers();
	auto dest_buffers = to->get_buffers();
	if(from_buffers.size() != dest_buffers.size()) {
		throw std::invalid_argument("The lattices need to have the same number of buffers.");
	}

	for(size_t i = 0; i < from_buffers.size(); ++i) {
		hardware::buffers::copyData(dest_buffers[i], from_buffers[i]);
	}
}

template<class T> void copyData(const T* to, const T* from)
{
	copyData(to, *from);
}

template <class Lattice, typename Basetype> void pseudo_randomize(const Lattice* to, int seed)
{
	unsigned elems = to->get_elements();
	std::vector<Basetype> host_vals(elems);
	Basetype * host_vals_p = host_vals.data();
	fill(host_vals_p, elems, seed);
	to->import(host_vals_p);
}

}

}

#endif /* _PHYSICS_LATTICES_UTIL_ */
