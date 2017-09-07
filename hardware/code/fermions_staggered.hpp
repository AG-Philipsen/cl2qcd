/** @file
 * Fermions_staggered OpenCL functionality
 *
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

#ifndef _HARDWARE_CODE_FERMIONS_STAGGERED_
#define _HARDWARE_CODE_FERMIONS_STAGGERED_

#include "opencl_module.hpp"

#include "../buffers/plain.hpp"
#include "../buffers/su3.hpp"
#include "../buffers/su3vec.hpp"
#include "../../host_functionality/host_use_timer.h"

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
class Fermions_staggered : public hardware::code::Opencl_Module {
public:
	friend hardware::Device;

	virtual ~Fermions_staggered();

	/*********************************************************************************************/
	/**************************  NON EVEN-ODD PRECONDITIONING METHODS  ***************************/
	/*********************************************************************************************/
	
	/**
	 * This function apply the standard staggered Dirac operator M = D_KS + m
	 * to a staggered field on the whole lattice.
	 *  @param in The input staggered field
	 *  @param out The output staggered field out=M*in
	 *  @param gf The gauge configuration
	 *  @param mass The quark mass
	 */
	void M_staggered_device(const hardware::buffers::Plain<su3vec> * in, const hardware::buffers::Plain<su3vec> * out, const hardware::buffers::SU3 * gf, hmc_float mass = ARG_DEF) const;

	/*********************************************************************************************/
	/****************************  EVEN-ODD PRECONDITIONING METHODS  *****************************/
	/*********************************************************************************************/
	
	/**
	 * This function apply the operator D_KS to a staggered field on half lattice.
	 * Actually, depending on the value of the variable evenodd, either the operator Doe
	 * or the operator Deo is selected and applied: if evenodd==EVEN then the staggered
	 * field in must be an ODD field and out=Deo*in is returned; if evenodd==ODD then 
	 * the staggered field in must be an EVEN field and out=Doe*in is returned. 
	 * @note Since there cannot be any check in the code, if one invoke this function with
	 *       evenodd equal to the parity of the field, there will be no error thrown.
	 *       By the programming point of you both an even field and an odd field are 
	 *       vectors with VOL4D/2 components.
	 * \par
	 * @note This function does not require the mass of the fermions as an input, because
	 *       in the staggered formulation the mass appear only in the diagonal part of the
	 *       Dirac operator.
	 *       
	 *  @param in The even/odd input staggered field
	 *  @param out The odd/even output staggered field out=D_KS*in
	 *  @param gf The gauge configuration
	 */
	void D_KS_eo_device(const hardware::buffers::SU3vec * in, const hardware::buffers::SU3vec * out, const hardware::buffers::SU3 * gf, int evenodd) const;
	
	////////////////////////////////////////////////////////////////////////////////////////
	/**
	 * Print the profiling information to a file.
	 *
	 * @param filename Name of file where data is appended.
	 */
	void virtual print_profiling(const std::string& filename, int number) const override;

	/**
	 * Return amount of bytes read and written by a specific kernel per call.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual size_t get_read_write_size(const std::string& in) const override;

	/**
	 * Return amount of Floating point operations performed by a specific kernel per call.
	 * @NOTE: this is meant to be the "netto" amount in order to be comparable.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual uint64_t get_flop_size(const std::string& in) const override;

	ClSourcePackage get_sources() const noexcept;

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
	 * @todo: the constructor must be public at the moment in order to be called from OpenClCode class.
	 * 	It may be made private again in the future!
	 */
public:
	Fermions_staggered(const hardware::code::OpenClKernelParametersInterface& kernelParams, const hardware::Device * device);

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

	////////////////////////////////////
	//fermionmatrix
	cl_kernel M_staggered;
	cl_kernel D_KS_eo;

	ClSourcePackage sources;
};

}

}

#endif // _HARDWARE_CODE_FERMIONS_STAGGERED_
