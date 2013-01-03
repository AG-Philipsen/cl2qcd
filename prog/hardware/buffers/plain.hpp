/** @file
 * Declaration of the hardware::buffers::Plain template
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
		template<typename T> class Plain : public Buffer {

		public:

			/**
			 * Create a buffer with the given number of Elements.
			 *
			 * \param elems The number of elements the buffer should contain
			 * \param device The device on which to create the buffer
			 */
			Plain(size_t elems, Device * device)
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
			 * Stores the whole buffer into the given pointer
			 *
			 * The value must not be used until the returned SynchronizationEvent
			 * returns true on is_finished().
			 */
			hardware::SynchronizationEvent dump_async(T * ptr) const
			{
				return Buffer::dump_async(ptr);
			}

			/**
			 * Get the size of the buffer in elements
			 */
			inline size_t get_elements() const noexcept
			{
				return elements;
			}

			/**
			 * Get data from another buffer.
			 * If the whole buffer should be copied use copyData(dest, src) instead!
			 *
			 * Will thorw an invalid_argument exception if the source buffer is of a different size.
			 */
			void copyDataBlock(const Buffer* orig, const size_t dest_offset, const size_t src_offset = 0, size_t elems = 0) const {
				if(!elems) {
					elems = this->elements;
				}
				Buffer::copyDataBlock(orig, dest_offset * sizeof(T), src_offset * sizeof(T), elems * sizeof(T));
			}

		private:
			const size_t elements;
		};
	}
}

#endif /* _HARDWARE_BUFFERS_SCALAR_BUFFER_ */
