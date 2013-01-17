/** @file
 * Functions for working with sources
 */

#ifndef _PHYSICS_SOURCES_
#define _PHYSICS_SOURCES_

#include "../hardware/system.hpp"
#include "prng.hpp"
#include "lattices/spinorfield.hpp"
#include "lattices/swappable_spinorfield.hpp"

namespace physics {

	/**
	 * Create sources as specified by the input parameters of the system.
	 */
	std::vector<lattices::Spinorfield *> create_sources(hardware::System& system, PRNG& prng);

	/**
	 * Create a set of spinorfields that can be swapped.
	 * Return normal spinorfield pointers for compatibility reasons.
	 */
	std::vector<lattices::Spinorfield *> create_swappable_sources(hardware::System& system, PRNG& prng);

	void set_point_source(const physics::lattices::Spinorfield *, int k, const meta::Inputparameters& params);
	void set_volume_source(const physics::lattices::Spinorfield *, PRNG& prng);
	void set_timeslice_source(const physics::lattices::Spinorfield *, PRNG& prng, int t);
	void set_zslice_source(const physics::lattices::Spinorfield *, PRNG& prng, int z);

}

#endif /* _PHYSICS_SOURCES_ */
