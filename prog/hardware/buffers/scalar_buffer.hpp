/** @file
 * Declaration of the hardware::buffers::ScalarBuffer template
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
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
		template<typename T> class ScalarBuffer : public Buffer {

		public:

			/**
			 * Create a buffer with the given number of Elements.
			 *
			 * \param elems The number of elements the buffer should contain
			 * \param device The device on which to create the buffer
			 */
			ScalarBuffer(size_t elems, Device * device)
				: Buffer(elems * sizeof(T), device), elements(elems) { };

			/**
			 * Loads as many elements as the buffer contains from the given pointer
			 * into the buffer.
			 */
			inline void load(const T* ptr) const
			{
				Buffer::load(ptr);
			}

			/**
			 * Stores the whole buffer into the given pointer
			 */
			inline void dump(T* ptr) const
			{
				Buffer::dump(ptr);
			}

			/**
			 * Get the size of the buffer in elements
			 */
			inline size_t get_elements() const noexcept
			{
				return elements;
			}

		private:
			const size_t elements;
		};
	}
}

#endif /* _HARDWARE_BUFFERS_SCALAR_BUFFER_ */
