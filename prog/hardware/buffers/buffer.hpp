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
			operator const cl_mem*() const noexcept;

			/**
			 * Also have a function to get the cl buffer,
			 * easier when having a pointer...
			 */
			const cl_mem* get_cl_buffer() const noexcept;
	
			/**
			 * Get the buffer size in bytes
			 */
			size_t get_bytes() const noexcept;

			/**
			 * Get the device this buffer is located on
			 */
			hardware::Device * get_device() const noexcept;

		private:
			/**
			 * The size of the buffer in bytes.
			 */
			const size_t bytes;

			/**
			 * The OpenCL buffer handle.
			 */
			const cl_mem cl_buffer;

			/**
			 * The OpenCL device the buffer is located on.
			 */
			Device * const device;

		protected:
			/**
			 * Utility function for creation of custom load functions.
			 *
			 * Loads as many bytes as the buffer contains from the given pointer
			 * into the buffer.
			 */
			void load(const void*) const;

			/**
			 * Utility function for creation of custom dump functions.
			 * 
			 * Stores the whole buffer into the given pointer
			 */
			void dump(void*) const;
		};
	}
}

#endif /* _HARDWARE_BUFFERS_BUFFER_ */
