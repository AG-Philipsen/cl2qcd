/** @file
 * Declaration of the hardware::buffers::SU3 class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#ifndef _HARDWARE_BUFFERS_SU3_
#define _HARDWARE_BUFFERS_SU3_

#include "buffer.hpp"
#include "../../types.h"

namespace hardware {
	namespace buffers {

		/**
		 * Check whether SU3 should be stored SOA style on this device
		 */
		size_t check_SU3_for_SOA(hardware::Device * device);

		/**
		 * Get the stride for an SU3 buffer of the given number of elements on the given device
		 */
		size_t get_SU3_buffer_stride(size_t elems, Device * device);

		/*
		 * A buffer storing SU3s on the device
		 */
		class SU3 : public Buffer {

		public:
			/**
			 * Allocate a buffer with the default number of
			 * elemets for this device.
			 *
			 * \param elems The size of the buffer in elements
			 * \param device The device to locate the buffer on
			 */
			SU3(size_t elems, Device * device);

			/*
			 * SU3 buffers cannot be copied
			 */
			SU3& operator=(const SU3&) = delete;
			SU3(const SU3&) = delete;
			SU3() = delete;

			/**
			 * Load data from the given pointer into the buffer.
			 *
			 * This only works for AoS-Buffers. If the buffer is a SoA buffer
			 * an std::logic_error will be thrown.
			 */
			void load(const Matrixsu3 *) const;

			/**
			 * Store data from the buffer into the given pointer.
			 *
			 * This only works for AoS-Buffers. If the buffer is a SoA buffer
			 * an std::logic_error will be thrown.
			 */
			void dump(Matrixsu3 *) const;

			/**
			 * Get the size of the buffer in elements
			 */
			size_t get_elements() const noexcept;

			/**
			 * Check whether this Buffer uses soa layout
			 */
			bool is_soa() const noexcept;

		private:

			/**
			 * The size of the buffer in bytes.
			 */
			const size_t elems;

			/**
			 * Whether the data is stored in a soa fashion
			 */
			const bool soa;
		};
	}
}

#endif /* _HARDWARE_BUFFERS_PRNG_BUFFER_ */

