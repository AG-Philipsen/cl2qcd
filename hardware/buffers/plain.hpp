/** @file
 * Declaration of the hardware::buffers::Plain template
 *
 * Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
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


#ifndef _HARDWARE_BUFFERS_SCALAR_BUFFER_
#define _HARDWARE_BUFFERS_SCALAR_BUFFER_

#include "buffer.hpp"

namespace hardware {
namespace buffers {

/**
 * A generic Buffer implementation for "scalar" types.
 *
 * This provides a full implementation of the Buffer interface
 * for scalar and scalar-like types.
 * A scalar-like type can be any type that should always be treated like a scalar, e.g. stored in a AoS fashion.
 * While this often should not be done for performance reasons it might be required for import/export code.
 */
template<typename T> class Plain : public Buffer {

public:

	/**
	 * Create a buffer with the given number of Elements.
	 *
	 * \param elems The number of elements the buffer should contain
	 * \param device The device on which to create the buffer
	 * \param place_on_host Request the buffer to remain on the host
	 *                      Allows to avoid memory limits if many buffers are required
	 */
	Plain(const size_t elems, const Device * device, bool place_on_host = false, cl_mem_flags additional_flags = 0)
		: Buffer(elems * sizeof(T), device, place_on_host, additional_flags), elements(elems) { };

	/**
	 * Loads as many elements as the buffer contains from the given pointer
	 * into the buffer.
	 *
	 * \param elems Allows to limit the number of elements loaded from the given pointer
	 * \param offset Allows to store the elements at the given offset into the buffer
	 */
	inline void load(const T* ptr, size_t elems = 0, size_t offset = 0) const {
		Buffer::load(ptr, elems * sizeof(T), offset * sizeof(T));
	}

	/**
	 * Stores the whole buffer into the given pointer
	 *
	 * \param elems Allows to limit the number of elements dumped to the given pointer
	 * \param offset Allows to read the elements at the given offset into the buffer
	 */
	inline void dump(T* ptr, size_t elems = 0, size_t offset = 0) const {
		Buffer::dump(ptr, elems * sizeof(T), offset * sizeof(T));
	}

	/**
	 * Loads the whole buffer into the given pointer
	 *
	 * The value must not be modified or deleted until the returned SynchronizationEvent
	 * returns true on is_finished().
	 */
	hardware::SynchronizationEvent load_async(const T * ptr) const {
		return Buffer::load_async(ptr);
	}

	/**
	 * Stores the whole buffer into the given pointer
	 *
	 * The value must not be used until the returned SynchronizationEvent
	 * returns true on is_finished().
	 */
	hardware::SynchronizationEvent dump_async(T * ptr) const {
		return Buffer::dump_async(ptr);
	}

	/**
	 * Get the size of the buffer in elements
	 */
	inline size_t get_elements() const noexcept {
		return elements;
	}

	/**
	 * Get data from another buffer.
	 * If the whole buffer should be copied use copyData(dest, src) instead!
	 */
	void copyDataBlock(const Plain<T>* orig, const size_t dest_offset, const size_t src_offset = 0, size_t elems = 0) const {
		if(!elems) {
			elems = std::min(this->elements, orig->elements);
		}
		Buffer::copyDataBlock(orig, dest_offset * sizeof(T), src_offset * sizeof(T), elems * sizeof(T));
	}

private:
	const size_t elements;
};
}
}

#endif /* _HARDWARE_BUFFERS_SCALAR_BUFFER_ */
