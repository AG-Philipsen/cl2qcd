/** @file
 * Declaration of the hardware::buffers::Buffer class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#ifndef _HARDWARE_BUFFERS_BUFFER_
#define _HARDWARE_BUFFERS_BUFFER_

#include "../device.hpp"

namespace hardware {
	/**
	 * This namespace contains all buffers, physically representing fields
	 * of specific data on the device.
	 */
	namespace buffers {
	
		/**
		 * A generic OpenCL buffer.
		 *
		 * This class should not really be used directly, it more works as a
		 * common RAII implementation for more specific buffer classes.
		 */
		class Buffer {
	
		public:
			/**
			 * Allocate a buffer of the given size on the given device.
			 *
			 * \param bytes The size of the buffer in bytes
			 * \param device The device to locate the buffer on
			 */
			Buffer(size_t bytes, Device * device);
	
			virtual ~Buffer();
	
			/*
			 * Buffers cannot be copied
			 */
			Buffer& operator=(const Buffer&) = delete;
			Buffer(const Buffer&) = delete;
			Buffer() = delete;
	
			/**
			 * Get a pointer to the opencl buffer object.
			 * This allows directly passing this class to
			 * clSetKernelArg
			 */
			operator const cl_mem*() const;
	
			/**
			 * Get the buffer size in bytes
			 */
			size_t get_bytes() const;

		private:
			/**
			 * The size of the buffer in bytes.
			 */
			const size_t bytes;

			/**
			 * The OpenCL buffer handle.
			 */
			const cl_mem cl_buffer;
		};
	}
}

#endif /* _HARDWARE_BUFFERS_BUFFER_ */
