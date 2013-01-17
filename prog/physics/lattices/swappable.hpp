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
