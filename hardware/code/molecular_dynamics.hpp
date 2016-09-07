/** @file
 * Molecular_Dynamics OpenCL functionality
 *
 * Copyright (c) 2012 Christopher Pinke <pinke@compeng.uni-frankfurt.de>
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.phys.uni-frankfurt.de>
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


#ifndef _HARDWARE_CODE_MOLECULAR_DYNAMICS_
#define _HARDWARE_CODE_MOLECULAR_DYNAMICS_

#include "opencl_module.hpp"
#include "../buffers/plain.hpp"
#include "../buffers/prng_buffer.hpp"
#include "../buffers/su3.hpp"
#include "../buffers/su3vec.hpp"
#include "../buffers/3x3.hpp"
#include "../buffers/spinor.hpp"
#include "../buffers/gaugemomentum.hpp"

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
class Molecular_Dynamics : public Opencl_Module {
public:
	friend hardware::Device;

	virtual ~Molecular_Dynamics();

	///////////////////////////////////////////////////
	//Methods on device

	void md_update_gaugefield_device(hmc_float eps) const;
	void md_update_gaugefield_device(const hardware::buffers::Gaugemomentum *, const hardware::buffers::SU3 *, hmc_float eps) const;
	void gauge_force_device(const hardware::buffers::SU3 * gf, const hardware::buffers::Gaugemomentum * out) const;
	void gauge_force_tlsym_device(const hardware::buffers::SU3 * gf, const hardware::buffers::Gaugemomentum * out) const;
	void fermion_force_device(const hardware::buffers::Plain<spinor> * Y, const hardware::buffers::Plain<spinor> * X, hmc_float kappa = ARG_DEF) const;
	void fermion_force_device(const hardware::buffers::Plain<spinor> * Y, const hardware::buffers::Plain<spinor> * X, const hardware::buffers::SU3 *, const hardware::buffers::Gaugemomentum *, hmc_float kappa = ARG_DEF) const;
	void fermion_force_eo_device(const hardware::buffers::Spinor * Y, const hardware::buffers::Spinor * X, int evenodd, hmc_float kappa = ARG_DEF) const;
	void fermion_force_eo_device(const hardware::buffers::Spinor * Y, const hardware::buffers::Spinor * X, const hardware::buffers::SU3 *, const hardware::buffers::Gaugemomentum *, int evenodd, hmc_float kappa = ARG_DEF) const;
	void stout_smeared_fermion_force_device(std::vector<const hardware::buffers::SU3 *>& gf_intermediate) const;
	///////////////////////////////////////////////////
	//Methods added exclusively for staggered fermions
	void fermion_staggered_partial_force_device(const hardware::buffers::SU3 * gf, const hardware::buffers::SU3vec * A, const hardware::buffers::SU3vec * B, const hardware::buffers::Gaugemomentum * out, int evenodd) const;

	/**
	 * Print the profiling information to a file.
	 *
	 * @param filename Name of file where data is appended.
	 */
	void virtual print_profiling(const std::string& filename, int number) const override;

protected:

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
	 * @todo: the constructor must be public at the moment in order to be called from OpenClCode class.
	 * 	It may be made private again in the future!
	 */
public:
	Molecular_Dynamics(const hardware::code::OpenClKernelParametersInterface& kernelParams, const hardware::Device * device);

private:
	/**
	 * Collect the kernels for OpenCL.
	 */
	void fill_kernels();
	/**
	 * Clear out the kernels,
	 */
	void clear_kernels();

	ClSourcePackage basic_molecular_dynamics_code;

	//kernels
	cl_kernel md_update_gaugefield;
	cl_kernel gauge_force;
	cl_kernel gauge_force_tlsym;
	cl_kernel fermion_force;
	cl_kernel fermion_force_eo_0;
	cl_kernel fermion_force_eo_1;
	cl_kernel fermion_force_eo_2;
	cl_kernel fermion_force_eo_3;
	cl_kernel stout_smear_fermion_force;
	
	//staggered kernels
	cl_kernel fermion_stagg_partial_force_eo;

	cl_kernel gauge_force_tlsym_1;
	cl_kernel gauge_force_tlsym_2;
	cl_kernel gauge_force_tlsym_3;
	cl_kernel gauge_force_tlsym_4;
	cl_kernel gauge_force_tlsym_5;
	cl_kernel gauge_force_tlsym_6;

	const hardware::buffers::Matrix3x3 * gauge_force_tlsym_tmp;
};

}

}

#endif // _HARDWARE_CODE_MOLECULAR_DYNAMICS_
