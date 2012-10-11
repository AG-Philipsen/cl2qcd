/** @file
 * Declaration of the hardware::buffers::Spinor class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#ifndef _HARDWARE_BUFFERS_SPINOR_
#define _HARDWARE_BUFFERS_SPINOR_

#include "buffer.hpp"
#include "../../types_fermions.h"

namespace hardware {
	namespace buffers {

		/**
		 * Check whether Spinor should be stored SOA style on this device
		 */
		size_t check_Spinor_for_SOA(hardware::Device * device);

		/**
		 * Get the stride for an Spinor buffer of the given number of elements on the given device
		 */
		size_t get_Spinor_buffer_stride(size_t elems, Device * device);

		/*
		 * A buffer storing Spinors on the device
		 */
		class Spinor : public Buffer {

		public:
			/**
			 * Allocate a buffer with the default number of
			 * elemets for this device.
			 *
			 * \param elems The size of the buffer in elements
			 * \param device The device to locate the buffer on
			 */
			Spinor(size_t elems, Device * device);

			/*
			 * Spinor buffers cannot be copied
			 */
			Spinor& operator=(const Spinor&) = delete;
			Spinor(const Spinor&) = delete;
			Spinor() = delete;

			/**
			 * Load data from the given pointer into the buffer.
			 *
			 * This only works for AoS-Buffers. If the buffer is a SoA buffer
			 * an std::logic_error will be thrown.
			 */
			void load(const spinor *) const;

			/**
			 * Store data from the buffer into the given pointer.
			 *
			 * This only works for AoS-Buffers. If the buffer is a SoA buffer
			 * an std::logic_error will be thrown.
			 */
			void dump(spinor *) const;

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

#endif /* _HARDWARE_BUFFERS_SPINOR_ */
