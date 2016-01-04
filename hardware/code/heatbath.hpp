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
#ifndef _HARDWARE_CODE_HEATBATH_
#define _HARDWARE_CODE_HEATBATH_

#include "opencl_module.hpp"
#include "../buffers/su3.hpp"
#include "../buffers/prng_buffer.hpp"

namespace hardware {

namespace code {

/**
 * An OpenCL device
 *
 * Adds heatbath to basic Opencl_Module class
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Heatbath : public hardware::code::Opencl_Module {
public:
	friend hardware::Device;

	virtual ~Heatbath();

	/**
	 * Perform one heatbath step.
	 */
	void run_heatbath(const hardware::buffers::SU3 * gaugefield, const hardware::buffers::PRNGBuffer * prng) const;

	/**
	 * Perform one overrelaxation step.
	 */
	void run_overrelax(const hardware::buffers::SU3 * gaugefield, const hardware::buffers::PRNGBuffer * prng) const;

	/**
	 * Add specific work_size determination for this child class
	 */
	virtual void get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const override;

	/**
	 * Print the profiling information to a file.
	 *
	 * @param filename Name of file where data is appended.
	 * @param parameters inputparameters
	 */
	void virtual print_profiling(const std::string& filename, int number) const override;

	/**
	 * @todo: the constructor must be public at the moment in order to be called from OpenClCode class.
	 * 	It may be made private again in the future!
	 */
public:
	Heatbath(const hardware::code::OpenClKernelParametersInterface& kernelParams, const hardware::Device * device);

private:
	/**
	 * Collect the kernels for OpenCL.
	 * Virtual method, allows to include more kernels in inherited classes.
	 */
	void fill_kernels();

	/**
	 * Clear out the kernels,
	 * Virtual method, allows to clear additional kernels in inherited classes.
	 */
	void clear_kernels();

	cl_kernel heatbath_odd;
	cl_kernel heatbath_even;
	cl_kernel overrelax_odd;
	cl_kernel overrelax_even;


	/**
	 * Return amount of bytes read and written by a specific kernel per call.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual size_t get_read_write_size(const std::string& in) const override;

	/**
	 * Return amount of Floating point operations performed by a specific kernel per call.
	 * NOTE: this is meant to be the "netto" amount in order to be comparable.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual uint64_t get_flop_size(const std::string& in) const override;


};

}

}

#endif // _HARDWARE_CODE_HEATBATH_
