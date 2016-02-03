/** @file
 * Declaration of the physics::lattices::SwappableSpinorfield class
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

#ifndef _PHYSICS_LATTICES_SWAPPABLESPINORFIELD_
#define _PHYSICS_LATTICES_SWAPPABLESPINORFIELD_

#include "spinorfield.hpp"
#include "swappable.hpp"
#include "../interfacesHandler.hpp"

/**
 * This namespace contains the lattices of the various kind,
 * that is storage of the lattice values as a whole.
 */
namespace physics {
namespace lattices {

/**
 * Representation of a gaugefield.
 */
class SwappableSpinorfield : public Spinorfield, public Swappable {
public:
	/**
	 * Construct a gaugefield based on the input-files of the system
	 */
	SwappableSpinorfield(const hardware::System&, const physics::lattices::SpinorfieldParametersInterface & spinorfieldParametersInterface, const bool place_on_host = false);

	/**
	 * Release resources
	 */
	virtual ~SwappableSpinorfield();

	/**
	 * Ensure the data is available for usage on the device.
	 *
	 * If the data is not available on the device simply fail over.
	 */
	void swap_in();

	/**
	 * Ensure the data does not take up space on the device.
	 *
	 * If the data has already been swapped_out simple fail over.
	 * Ensures that a call to get_buffers() will not return any buffers.
	 */
	void swap_out();

private:
	// TODO list would probably be better
	std::vector<spinor*> swap;
};

/**
 * Create n spinorfields that can be swapped.
 * For combatibiltiy reasons return normal spinorfield pointers.
 *
 * \param n The number of spinorfields to create
 */
std::vector<Spinorfield *> create_swappable_spinorfields(const hardware::System& system, const size_t n, physics::InterfacesHandler & interfacesHandler, const bool place_on_host = false);

/**
 * Swap out the given spinorfields
 *
 * \throws Print_Error_Message if not all Spinorfields are swappable
 */
void swap_out(const std::vector<Spinorfield *>& fields);

/**
 * Swap in the given spinorfields
 *
 * If a spinorfield is not swappable this is not an error.
 */
void swap_in(const std::vector<Spinorfield *>& fields);

}
}

#endif /*_PHYSICS_LATTICES_SWAPPABLESPINORFIELD_ */
