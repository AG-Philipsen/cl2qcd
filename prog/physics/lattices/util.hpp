/** @file
 * Utility functions for physics::lattices
 */

#ifndef _PHYSICS_LATTICES_UTIL_
#define _PHYSICS_LATTICES_UTIL_

#include <stdexcept>
#include "../../hardware/buffers/buffer.hpp"
#include "../../meta/type_ops.hpp"

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
template <class T> void pseudo_randomize(int seed);



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
	auto buffers = to->get_buffers();
	if(buffers.size() != 1) {
		throw Print_Error_Message("Pseudo-randomization of multi-buffer lattices is not yet implemented.");
	}
	auto buffer = buffers[0];
	size_t elems = buffer->get_elements();
	std::vector<Basetype> host_vals(elems);
	Basetype * host_vals_p = host_vals.data();
	fill(host_vals_p, elems, seed);
	buffer->load(host_vals_p);
}

}

}

#endif /* _PHYSICS_LATTICES_UTIL_ */
