/** @file
 * Utility functions for physics::lattices
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



/*
 * TEMPLATE IMPLEMENTATIONS
 */
template<class T> void copyData(const T* to, const T* from)
{
	auto from_buffers = from->get_buffers();
	auto dest_buffers = to->get_buffers();
	if(from_buffers.size() != dest_buffers.size()) {
		throw std::invalid_argument("The lattices need to have the same number of buffers.");
	}

	for(size_t i = 0; i < from_buffers.size(); ++i) {
		hardware::buffers::copyData(dest_buffers[i], from_buffers[i]);
	}
}

}

}

#endif /* _PHYSICS_LATTICES_UTIL_ */
