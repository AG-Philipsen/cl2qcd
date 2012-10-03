/** @file
 * Declaration of the hardware::buffers::Buffer class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#ifndef _HARDWARE_BUFFERS_PRN_GBUFFER_
#define _HARDWARE_BUFFERS_PRN_GBUFFER_

#include "buffer.hpp"

namespace hardware {
	namespace buffers {

		/**
		 * Get the default size PRNG buffer on the given device.
		 *
		 * \param device The device the buffer is to be used with
		 * \return The size of the buffer in element
		 */
		size_t get_prng_buffer_size(const Device * device);

		/**
		 * A PRNG OpenCL buffer.
		 */
		class PRNGBuffer : public Buffer {

		public:
#ifdef USE_PRNG_NR3
			typedef nr3_state_dev prng_state_t;
#elif defined(USE_PRNG_RANLUX)
			typedef cl_float4 prng_state_t[7];
			static_assert(sizeof(prng_state_t) == 7 * sizeof(cl_float4), "PRNG state type mockup is of wrong size.");
#else // USE_PRNG_XXX
#error No implemented PRNG selected
#endif // USE_PRNG_XXX

			/**
			 * Allocate a buffer with the default number of
			 * elemets for this device.
			 *
			 * \param elems The size of the buffer in elements
			 * \param device The device to locate the buffer on
			 */
			PRNGBuffer(Device * device);

			/*
			 * PRNGBuffers cannot be copied
			 */
			PRNGBuffer& operator=(const PRNGBuffer&) = delete;
			PRNGBuffer(const PRNGBuffer&) = delete;
			PRNGBuffer() = delete;

			/**
			 * Get the buffer size in bytes
			 */
			size_t get_elements() const noexcept;

			void load(const prng_state_t *) const;

			void dump(prng_state_t *) const;

		private:
			/**
			 * Allocate a buffer of the given size in elements.
			 *
			 * \param elems The size of the buffer in elements
			 * \param device The device to locate the buffer on
			 */
			PRNGBuffer(size_t elems, Device * device);

			/**
			 * The size of the buffer in bytes.
			 */
			const size_t elems;
		};
	}
}

#endif /* _HARDWARE_BUFFERS_PRNG_BUFFER_ */

