/** @file
 * Basic OpenCL functionality
 *
 * Copyright 2012, 2013, 2016 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra,
 * Tim Breitenfelder
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

#ifndef _HARDWARE_CODE_CORRELATOR_STAGGERED_
#define _HARDWARE_CODE_CORRELATOR_STAGGERED_

#include "opencl_module.hpp"

#include "../buffers/plain.hpp"
#include "../buffers/su3vec.hpp"
#include "../buffers/prng_buffer.hpp"
#include "../../common_header_files/types_fermions.h"

namespace hardware {

namespace code {

/**
 * An OpenCL device
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Correlator_staggered : public Opencl_Module {
public:
	friend hardware::Device;

	virtual ~Correlator_staggered();

	void create_point_source_stagg_eoprec_device(const hardware::buffers::SU3vec * inout, int i, int spacepos, int timepos) const;

	void create_volume_source_stagg_eoprec_device(const hardware::buffers::SU3vec * inout, const hardware::buffers::PRNGBuffer * prng) const;

    /**
     * Calculate the correlator on the device.
     * TODO: In future, if different correlators are added, think whether to do as in Wilson, overloading, or in a different way!
     */
    void pseudoScalarCorrelator(const hardware::buffers::Plain<hmc_float> * correlator, const hardware::buffers::SU3vec * invertedSourcesEven, const hardware::buffers::SU3vec * invertedSourcesOdd) const;

	/**
	 * Print the profiling information to a file.
	 *
	 * @param filename Name of file where data is appended.
	 * @param number task-id
	 */
	void virtual print_profiling(const std::string& filename, int number) const override;

protected:
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

	/**
	 * comutes work-sizes for a kernel
	 * @todo autotune
	 * @param ls local-work-size
	 * @param gs global-work-size
	 * @param num_groups number of work groups
	 * @param name name of the kernel for possible autotune-usage, not yet used!!
	 */
	virtual void get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const override;

	/**
	 * @todo: the constructor must be public at the moment in order to be called from OpenClCode class.
	 * 	It may be made private again in the future!
	 */
public:
	Correlator_staggered(const hardware::code::OpenClKernelParametersInterface& kernelParams , const hardware::Device * device);

private:
	/**
	 * Collect the kernels for OpenCL.
	 */
	void fill_kernels();
	/**
	 * Clear out the kernels,
	 */
	void clear_kernels();

	////////////////////////////////////
	//kernels
	///////////////////////////////////
	//Types of source

	cl_kernel create_point_source_stagg_eoprec;
	cl_kernel create_volume_source_stagg_eoprec;

	//Observables
	//scalar correlators
	cl_kernel correlator_staggered_ps;

	ClSourcePackage basic_correlator_code;
};

}

}

#endif // _HARDWARE_CODE_CORRELATOR_STAGGERED_
