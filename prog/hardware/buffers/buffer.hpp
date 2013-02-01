/** @file
 * Declaration of the hardware::buffers::Buffer class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#ifndef _HARDWARE_BUFFERS_BUFFER_
#define _HARDWARE_BUFFERS_BUFFER_

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <stdexcept>
#include "../synchronization_event.hpp"

namespace hardware {

	/* forward decleration of device to avoid cycles */
	class Device;

	/**
	 * This namespace contains all buffers, physically representing fields
	 * of specific data on the device.
	 */
	namespace buffers {

		/**
		 * Copy buffer data
		 *
		 * A utility function to copy buffer contents.
		 * This will throw an exception if the two buffers are not of equal size.
		 *
		 * \param dest The buffer to copy to
		 * \param src  The buffer to copy from
		 */
		template<class T> inline void copyData(const T* dest, const T* orig);

		/**
		 * A generic OpenCL buffer.
		 *
		 * This class should not really be used directly, it more works as a
		 * common RAII implementation for more specific buffer classes.
		 */
		class Buffer {

			template<class T> friend void copyData(const T* dest, const T* orig);
	
		public:
			/**
			 * Allocate a buffer of the given size on the given device.
			 *
			 * \param bytes The size of the buffer in bytes
			 * \param device The device to locate the buffer on
			 * \param place_on_host Request the buffer to remain on the host
			 *                      Allows to avoid memory limits if many buffers are required
			 */
			Buffer(size_t bytes, Device * device, bool place_on_host = false);
	
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

			/**
			 * Set all bytes of this buffer to zero.
			 */
			void clear() const;

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

			/**
			 * Utility function to get the data from another buffer. Should only be used using
			 * the copyData wrapper template to ensure proper type checking.
			 *
			 * Will thorw an invalid_argument exception if the source buffer is of a different size.
			 */
			void copyData(const Buffer* orig) const;

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

			/**
			 * Utility function for creation of custom dump functions.
			 *
			 * Stores the whole buffer into the given pointer
			 *
			 * The value must not be used until the returned SynchronizationEvent
			 * returns true on is_finished().
			 */
			hardware::SynchronizationEvent dump_async(void * array) const;

			/**
			 * Utility function to get the data from another buffer. Should only be used using
			 * Should be implemented by children using element instead of bytes sizes
			 *
			 * Will thorw an invalid_argument exception if the source buffer is of a different size.
			 */
			void copyDataBlock(const Buffer* orig, const size_t dest_offset, const size_t src_offset, const size_t bytes) const;

			friend std::string md5(const Buffer* buf);
		};

		template<class T> inline void copyData(const T* dest, const T* orig)
		{
			dest->copyData(orig);
		}

		/**
		 * Get the md5 sum of the given buffer
		 */
		std::string md5(const Buffer* buf);
	}
}

/**
 * Disallow to use Buffers as kernel arguments.
 *
 * The contained cl_mem must be passed to the kernel.
 */
cl_int clSetKernelArg(cl_kernel, cl_uint, size_t, const hardware::buffers::Buffer *) = delete;
// the following would be the nicer solution, as it avoids problems due to ambiguities, but it causes issues with clang on OSX and boost
//#include <type_traits>
//template<typename T> inline cl_int clSetKernelArg(cl_kernel kernel, cl_uint idx, size_t bytes, const T * arg)
//{
//	// Check that we don't pass managment objects
//	static_assert(!std::is_convertible<T, const cl_mem *>::value, "You should pass a pointer to cl_mem, not the managment object itself to the kernel.");
//	// Invoke the original OpenCL function
//	clSetKernelArg(kernel, idx, bytes, static_cast<const void *>(arg));
//}

#endif /* _HARDWARE_BUFFERS_BUFFER_ */
