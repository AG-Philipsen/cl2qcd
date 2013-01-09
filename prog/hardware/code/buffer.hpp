/** @file
 * Buffer kernels
 */
#ifndef _HARDWARE_CODE_BUFFER_
#define _HARDWARE_CODE_BUFFER_

#include "opencl_module.hpp"

#include "../buffers/buffer.hpp"

namespace hardware {

namespace code {

/**
 * An OpenCL device
 *
 * Adds random numbers to basic Opencl_Module class
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Buffer : public Opencl_Module {
public:
	friend hardware::Device;

	virtual ~Buffer();

	/**
	 * Copies a 16 byte buffer.
	 *
	 * @warn Use hardware::buffers::copyData!
	 */
	void copy_16_bytes(const hardware::buffers::Buffer * dest, const hardware::buffers::Buffer * orig) const;

	/**
	 * Clears the given buffer
	 *
	 * \dest the buffer to set to zero
	 */
	void clear(const hardware::buffers::Buffer * dest) const;

protected:
	/**
	 * Return amount of Floating point operations performed by a specific kernel per call.
	 * NOTE: this is meant to be the "netto" amount in order to be comparable.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual uint64_t get_flop_size(const std::string&) const { return 0; };

	/**
	 * Return amount of bytes read and written by a specific kernel per call.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual size_t get_read_write_size(const std::string&) const { return 0; };

private:

	/**
	 * @param[in] params points to an instance of inputparameters
	 */
	Buffer(const meta::Inputparameters& params, hardware::Device * device);

	cl_kernel _copy_16_bytes;
	cl_kernel _clear_bytes;
	cl_kernel _clear_float4;
};

}

}

#endif // _HARDWARE_CODE_BUFFER_

