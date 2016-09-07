/** @file
 * Declaration of the hardware::buffers::Buffer class
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

#ifndef _HARDWARE_BUFFERS_BUFFER_
#define _HARDWARE_BUFFERS_BUFFER_

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <stdexcept>
#include "../synchronization_event.hpp"
#include <memory>

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
 * Handle to a buffer mapped to CPU memory.
 */
class MappedBufferHandle;

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
	 * \param any additional flags required for this buffer
	 */
	Buffer(const size_t bytes, const Device * device, const bool place_on_host = false, const cl_mem_flags extra_flags = 0);

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
	const hardware::Device * get_device() const noexcept;

	/**
	 * Set all bytes of this buffer to zero.
	 */
	void clear() const;

	/**
	 * Map the buffer into host address space.
	 */
	std::unique_ptr<MappedBufferHandle> map(cl_map_flags flags = CL_MAP_READ | CL_MAP_WRITE) const;

#ifdef CL_VERSION_1_2
	/**
	 * Migrate the bufer to a different device.
	 */
	void migrate(hardware::Device * device, const std::vector<hardware::SynchronizationEvent>& events, cl_mem_migration_flags flags = 0);
#endif

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
	const Device * device;

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
	 * Load the given number of bytes from the given buffer.
	 *
	 * \param src, The memory area to read from
	 * \param bytes The number of bytes to rea, will be the buffer size if 0d
	 * \param offset Offset into the buffer, at which to store the bytes
	 * \throws out_of_range if bytes + offset exceeds the buffer size
	 */
	void load(const void* src, size_t bytes = 0, const size_t offset = 0) const;

	/**
	 * Utility function for creation of custom dump functions.
	 *
	 * Store the given number of bytes into the given pointer.
	 *
	 * \param dest The memory area to store to
	 * \param bytes The number of bytes to write, will be the buffer size if 0
	 * \param offset Offset into the buffer, at which to store the bytes
	 * \throws out_of_range if bytes + offset exceeds the buffer size
	 */
	void dump(void* dest, size_t bytes = 0, const size_t offset = 0) const;

	void load_rect(const void* src, const size_t *buffer_origin, const size_t *host_origin, const size_t *region, size_t buffer_row_pitch, size_t buffer_slice_pitch, size_t host_row_pitch, size_t host_slice_pitch) const;

	void dump_rect(void* dest, const size_t *buffer_origin, const size_t *host_origin, const size_t *region, size_t buffer_row_pitch, size_t buffer_slice_pitch, size_t host_row_pitch, size_t host_slice_pitch) const;

	hardware::SynchronizationEvent load_rectAsync(const void* src, const size_t *buffer_origin, const size_t *host_origin, const size_t *region, size_t buffer_row_pitch, size_t buffer_slice_pitch, size_t host_row_pitch, size_t host_slice_pitch, const hardware::SynchronizationEvent& event) const;

	hardware::SynchronizationEvent dump_rectAsync(void* dest, const size_t *buffer_origin, const size_t *host_origin, const size_t *region, size_t buffer_row_pitch, size_t buffer_slice_pitch, size_t host_row_pitch, size_t host_slice_pitch, const hardware::SynchronizationEvent& event) const;

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
	 * Utility function for creation of custom dump functions.
	 *
	 * Stores the whole buffer into the given pointer
	 *
	 * The value must not be modified or deleted until the returned SynchronizationEvent
	 * returns true on is_finished().
	 */
	hardware::SynchronizationEvent load_async(const void * array) const;

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

/**
 * Utility function to get the data from another buffer. Should only be used using
 * Should be implemented by children using element instead of bytes sizes
 *
 * Will thorw an invalid_argument exception if the source buffer is of a different size.
 */
hardware::SynchronizationEvent copyDataRect(const hardware::Device* device, const hardware::buffers::Buffer* dest, const hardware::buffers::Buffer* orig, const size_t *dest_origin, const size_t *src_origin, const size_t *region, size_t dest_row_pitch, size_t dest_slice_pitch, size_t src_row_pitch, size_t src_slice_pitch, const std::vector<hardware::SynchronizationEvent>& events);

class MappedBufferHandle {
	friend Buffer;
	MappedBufferHandle(cl_mem buf, cl_command_queue queue, void * mapped_ptr, hardware::SynchronizationEvent map_event);
public:
	MappedBufferHandle(MappedBufferHandle const &) = delete;
	~MappedBufferHandle();
	void * get_mapped_ptr() const;
	hardware::SynchronizationEvent get_map_event() const;
private:
	cl_mem const buf;
	cl_command_queue const queue;
	void * const mapped_ptr;
	hardware::SynchronizationEvent const map_event;
};

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
