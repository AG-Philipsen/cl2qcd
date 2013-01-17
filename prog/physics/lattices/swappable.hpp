/** @file
 * Declaration of the physics::lattices::Swappable interface
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

}
}

#endif /* _PHYSICS_LATTICES_SWAPPABLE_ */
