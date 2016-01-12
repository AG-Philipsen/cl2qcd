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

#ifndef _HARDWARE_CODE_FERMIONS_
#define _HARDWARE_CODE_FERMIONS_

#include "opencl_module.hpp"

#include "../buffers/plain.hpp"
#include "../buffers/su3.hpp"
#include "../buffers/spinor.hpp"
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
class Fermions : public hardware::code::Opencl_Module {
public:
	friend hardware::Device;

	virtual ~Fermions();


	//    fermionmatrix operations
	//    non-eo
	//        explicit
	void M_wilson_device(const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF) const;
	void M_tm_plus_device(const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const;
	void M_tm_minus_device(const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const;
	void gamma5_device(const hardware::buffers::Plain<spinor> * inout) const;
	//    eo
	//        explicit
	void gamma5_eo_device(const hardware::buffers::Spinor * inout) const;
	void M_tm_inverse_sitediagonal_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, hmc_float mubar = ARG_DEF) const;
	void M_tm_sitediagonal_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, hmc_float mubar = ARG_DEF) const;
	void M_tm_inverse_sitediagonal_minus_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, hmc_float mubar = ARG_DEF) const;
	void M_tm_sitediagonal_minus_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, hmc_float mubar = ARG_DEF) const;
	/**
	 * Perform dslash_eo on the whole buffer.
	 */
	void dslash_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, int evenodd, hmc_float kappa = ARG_DEF) const;
	/**
	 * Perform dslash_eo only on the boundary sites that require the halo.
	 */
	void dslash_eo_boundary(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, int evenodd, hmc_float kappa = ARG_DEF) const;
	/**
	 * Perform dslash_eo only on the inner sites which do not require the halo.
	 */
	void dslash_eo_inner(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, int evenodd, hmc_float kappa = ARG_DEF) const;
	//        merged
//	void Aee_AND_gamma5_eo(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
//	void Aee_minus_AND_gamma5_eo(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void dslash_AND_M_tm_inverse_sitediagonal_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, int evenodd, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const;
	void dslash_AND_M_tm_inverse_sitediagonal_minus_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, const hardware::buffers::SU3 * gf, int evenodd, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF) const;
	void M_tm_sitediagonal_AND_gamma5_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, hmc_float mubar = ARG_DEF) const;
	void M_tm_sitediagonal_minus_AND_gamma5_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::Spinor * out, hmc_float mubar = ARG_DEF) const;
	void saxpy_AND_gamma5_eo_device(const hardware::buffers::Spinor * x, const hardware::buffers::Spinor * y, const hmc_complex alpha, const hardware::buffers::Spinor * out) const;

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
	 * NOTE: this is meant to be the "netto" amount in order to be comparable.
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
	Fermions(const hardware::code::OpenClKernelParametersInterface& kernelParams, const hardware::Device * device);

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
	//kernels, sorted roughly by groups
	//fermionmatrix
	cl_kernel M_wilson;
	cl_kernel gamma5;
	cl_kernel M_tm_plus;
	cl_kernel M_tm_minus;
	cl_kernel gamma5_eo;
	cl_kernel M_tm_sitediagonal;
	cl_kernel M_tm_inverse_sitediagonal;
	cl_kernel M_tm_sitediagonal_minus;
	cl_kernel M_tm_inverse_sitediagonal_minus;
	cl_kernel dslash_eo;
	cl_kernel _dslash_eo_boundary;
	cl_kernel _dslash_eo_inner;
	cl_kernel dslash_AND_M_tm_inverse_sitediagonal_eo;
	cl_kernel dslash_AND_M_tm_inverse_sitediagonal_minus_eo;
	cl_kernel M_tm_sitediagonal_AND_gamma5_eo;
	cl_kernel M_tm_sitediagonal_minus_AND_gamma5_eo;
	cl_kernel saxpy_AND_gamma5_eo;

	ClSourcePackage sources;
};

}

}

#endif // _HARDWARE_CODE_FERMIONS_
