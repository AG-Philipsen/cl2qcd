/** @file
 * Buffer kernels
 *
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
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
#ifndef _HARDWARE_CODE_BUFFER_
#define _HARDWARE_CODE_BUFFER_

#include "opencl_module.hpp"

#include "../buffers/buffer.hpp"
#include "gaugefield.hpp"

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
	virtual uint64_t get_flop_size(const std::string&) const {
		return 0;
	};

	/**
	 * Return amount of bytes read and written by a specific kernel per call.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual size_t get_read_write_size(const std::string&) const {
		return 0;
	};


	/**
	 * @todo: the constructor must be public at the moment in order to be called from OpenClCode class.
	 * 	It may be made private again in the future!
	 */
public:
	Buffer(const hardware::code::OpenClKernelParametersInterface& kernelParams, const hardware::Device * device);

private:
	cl_kernel _copy_16_bytes;
	cl_kernel _clear_bytes;
	cl_kernel _clear_float4;
};

}

}

#endif // _HARDWARE_CODE_BUFFER_

