/** @file
 * Heatbath for OpenCL
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
#ifndef _HARDWARE_CODE_PRNG_
#define _HARDWARE_CODE_PRNG_

#include "opencl_module.hpp"
#include "../../host_functionality/logger.hpp"
#include "../device.hpp"
#include <sstream>
#include "gaugefield.hpp"
#include "../buffers/prng_buffer.hpp"

namespace hardware {

namespace code {

/**
 * An OpenCL device
 *
 * Adds random numbers to basic Opencl_Module class
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Prng : public Opencl_Module {
public:
	friend hardware::Device;

	virtual ~Prng();

	ClSourcePackage get_sources() const noexcept;

#ifdef USE_PRNG_RANLUX
	/**
	 * Initialize the state of the PRNG with the given seed.
	 */
	void initialize(const hardware::buffers::PRNGBuffer * buffer, cl_uint seed) const;
#endif /* USE_PRNG_RANLUX */

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
	Prng(const hardware::code::OpenClKernelParametersInterface& kernelParams, const hardware::Device * device);

private:
	/**
	 * A set of sources required to use the PRNG.
	 */
	ClSourcePackage prng_code;

#ifdef USE_PRNG_RANLUX
	cl_kernel init_kernel;
#endif // USE_PRNG_???
};

}

}

#endif // _HARDWARE_CODE_PRNG_
