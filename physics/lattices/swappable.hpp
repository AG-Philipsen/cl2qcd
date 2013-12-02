/** @file
 * Declaration of the physics::lattices::Swappable interface
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

#ifndef _PHYSICS_LATTICES_SWAPPABLE_
#define _PHYSICS_LATTICES_SWAPPABLE_

namespace physics {
namespace lattices {

/**
 * Representation of a gaugefield.
 */
struct Swappable {

	/**
	 * Ensure the data is available for usage on the device.
	 *
	 * If the data is not available on the device simply fail over.
	 */
	virtual void swap_in() = 0;

	/**
	 * Ensure the data does not take up space on the device.
	 *
	 * If the data has already been swapped_out simple fail over.
	 * Ensures that a call to get_buffers() will not return any buffers.
	 */
	virtual void swap_out() = 0;
};

/**
 * Swap in the given field if it is a Swapable.
 */
template <class C> void try_swap_in(C* field);
/**
 * Swap out the given field if it is a Swapable.
 */
template <class C> void try_swap_out(C* field);
}
}

/*
 * Template implementation
 */
template <class C> void physics::lattices::try_swap_in(C* field)
{
	auto swappable = dynamic_cast<Swappable*>(field);
	if(swappable) {
		swappable->swap_in();
	}
}
template <class C> void physics::lattices::try_swap_out(C* field)
{
	auto swappable = dynamic_cast<Swappable*>(field);
	if(swappable) {
		swappable->swap_out();
	}
}

#endif /* _PHYSICS_LATTICES_SWAPPABLE_ */
