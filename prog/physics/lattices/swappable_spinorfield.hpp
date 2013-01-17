/** @file
 * Declaration of the physics::lattices::SwappableSpinorfield class
 */

#ifndef _PHYSICS_LATTICES_SWAPPABLESPINORFIELD_
#define _PHYSICS_LATTICES_SWAPPABLESPINORFIELD_

#include "spinorfield.hpp"
#include "swappable.hpp"

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
	SwappableSpinorfield(const hardware::System&, const bool place_on_host = false);

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

}
}

#endif /*_PHYSICS_LATTICES_SWAPPABLESPINORFIELD_ */

