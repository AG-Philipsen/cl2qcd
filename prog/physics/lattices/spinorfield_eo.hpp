/** @file
 * Declaration of the physics::lattices::Spinorfield_eo class
 */

#ifndef _PHYSICS_LATTICES_SPINORFIELD_EO_
#define _PHYSICS_LATTICES_SPINORFIELD_EO_

#include "../../hardware/system.hpp"
#include "../../hardware/buffers/spinor.hpp"
#include "../prng.hpp"

/**
 * This namespace contains the lattices of the various kind,
 * that is storage of the lattice values as a whole.
 */
namespace physics {
namespace lattices {

/**
 * Representation of a gaugefield.
 */
class Spinorfield_eo {

public:
	/**
	 * Construct a gaugefield based on the input-files of the system
	 */
	Spinorfield_eo(hardware::System&);

	/**
	 * Release resources
	 */
	~Spinorfield_eo();

	/*
	 * Spinorfield_eos cannot be copied
	 */
	Spinorfield_eo& operator=(const Spinorfield_eo&) = delete;
	Spinorfield_eo(const Spinorfield_eo&) = delete;
	Spinorfield_eo() = delete;

	/**
	 * Get the buffers containing the gaugefield state on the devices.
	 */
	const std::vector<const hardware::buffers::Spinor *> get_buffers() const noexcept;

private:
	hardware::System const& system;
	const std::vector<const hardware::buffers::Spinor *> buffers;
};
}
}

#endif /*_PHYSICS_LATTICES_SPINORFIELD_EO_ */
