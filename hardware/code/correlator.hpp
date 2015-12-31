/** @file
 * Basic OpenCL functionality
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

#ifndef _HARDWARE_CODE_CORRELATOR_
#define _HARDWARE_CODE_CORRELATOR_

#include "opencl_module.hpp"

#include "../buffers/plain.hpp"
#include "../buffers/prng_buffer.hpp"
#include "../../common_header_files/types_fermions.h"

namespace hardware {

namespace code {

/**
 * An OpenCL device
 *
 * This class wraps all operations on a device. Operations are always specific, e.g. each kernel and copy operation
 * has it's own wrapper function.
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Correlator : public Opencl_Module {
public:
	friend hardware::Device;

	virtual ~Correlator();

	void create_point_source_device(const hardware::buffers::Plain<spinor> * inout, int i, int spacepos, int timepos) const;

	void create_volume_source_device(const hardware::buffers::Plain<spinor> * inout, const hardware::buffers::PRNGBuffer * prng) const;

	void create_timeslice_source_device(const hardware::buffers::Plain<spinor> * inout, const hardware::buffers::PRNGBuffer * prng, const int timeslice) const;

	void create_zslice_source_device(const hardware::buffers::Plain<spinor> * inout, const hardware::buffers::PRNGBuffer * prng, const int zslice) const;

	/**
	 * Calculate the correlator on the device.
	 * This function is overloaded depending on whether one needs the source for the calculation or not.
	 */
	void correlator(const cl_kernel correlator_kernel, const hardware::buffers::Plain<hmc_float> * correlator, const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * source = nullptr) const;
	void correlator(const cl_kernel correlator_kernel, const hardware::buffers::Plain<hmc_float> * correlator, const hardware::buffers::Plain<spinor> * in1, const hardware::buffers::Plain<spinor> * in2, const hardware::buffers::Plain<spinor> * in3, const hardware::buffers::Plain<spinor> * in4) const;
	void correlator(const cl_kernel correlator_kernel, const hardware::buffers::Plain<hmc_float> * correlator, const hardware::buffers::Plain<spinor> * in1, const hardware::buffers::Plain<spinor> * source1, const hardware::buffers::Plain<spinor> * in2, const hardware::buffers::Plain<spinor> * source2, const hardware::buffers::Plain<spinor> * in3, const hardware::buffers::Plain<spinor> * source3, const hardware::buffers::Plain<spinor> * in4, const hardware::buffers::Plain<spinor> * source4) const;

	/**
	 * Get kernel for correlator indicated by which
	 * @param[in] which string that identifies the correlator (ps or sc, vx, vy, vz, ax, ay, az)
	 * @return correlator_kernel
	 */
	cl_kernel get_correlator_kernel(std::string which) const;

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

public:
	Correlator(const hardware::code::OpenClKernelParametersInterface& kernelParams , const hardware::Device * device);

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

	cl_kernel create_point_source;
	cl_kernel create_volume_source;
	cl_kernel create_timeslice_source;
	cl_kernel create_zslice_source;

	//Observables
	//scalar correlators
	cl_kernel correlator_ps;
	cl_kernel correlator_sc;
	//vector correlators, physical basis: Gamma_4 * Gamma_{x,y,z}
	cl_kernel correlator_vx;
	cl_kernel correlator_vy;
	cl_kernel correlator_vz;
	//axial vector correlators, physical basis: Gamma_5 * Gamma_4 * Gamma_{x,y,z}
	cl_kernel correlator_ax;
	cl_kernel correlator_ay;
	cl_kernel correlator_az;
	//axial-vector pseudoscalar correlator
	cl_kernel correlator_avps;
	//chiral condensate
	cl_kernel pbp_std;
	cl_kernel pbp_tm_one_end;

	ClSourcePackage basic_correlator_code;
};

}

}

#endif // _HARDWARE_CODE_CORRELATOR_
