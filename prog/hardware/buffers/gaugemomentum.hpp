/** @file
 * Declaration of the hardware::buffers::Gaugemomentum class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#ifndef _HARDWARE_BUFFERS_GAUGEMOMENTUM_
#define _HARDWARE_BUFFERS_GAUGEMOMENTUM_

#include "buffer.hpp"
#include "../../types.h"

namespace hardware {
	namespace buffers {

		/**
		 * Check whether Gaugemomentum should be stored SOA style on this device
		 */
		size_t check_Gaugemomentum_for_SOA(hardware::Device * device);

		/**
		 * Get the stride for an Gaugemomentum buffer of the given number of elements on the given device
		 */
		size_t get_Gaugemomentum_buffer_stride(size_t elems, Device * device);

		/*
		 * A buffer storing Gaugemomentums on the device
		 */
		class Gaugemomentum : public Buffer {

		public:
			/**
			 * Allocate a buffer with the default number of
			 * elemets for this device.
			 *
			 * \param elems The size of the buffer in elements
			 * \param device The device to locate the buffer on
			 */
			Gaugemomentum(size_t elems, Device * device);

			/*
			 * Gaugemomentum buffers cannot be copied
			 */
			Gaugemomentum& operator=(const Gaugemomentum&) = delete;
			Gaugemomentum(const Gaugemomentum&) = delete;
			Gaugemomentum() = delete;

			/**
			 * Load data from the given pointer into the buffer.
			 *
			 * This only works for AoS-Buffers. If the buffer is a SoA buffer
			 * an std::logic_error will be thrown.
			 */
			void load(const ae *) const;

			/**
			 * Store data from the buffer into the given pointer.
			 *
			 * This only works for AoS-Buffers. If the buffer is a SoA buffer
			 * an std::logic_error will be thrown.
			 */
			void dump(ae *) const;

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

#endif /* _HARDWARE_BUFFERS_GAUGEMOMENTUM_ */
