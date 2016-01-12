/** @file
 * Gaugemomentum OpenCL functionality
 *
 * Copyright (c) 2012 Christopher Pinke <pinke@compeng.uni-frankfurt.de>
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.phys.uni-frankfurt.de>
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


#ifndef _HARDWARE_CODE_GAUGEMOMENTUM_
#define _HARDWARE_CODE_GAUGEMOMENTUM_

#include "opencl_module.hpp"
#include "../buffers/plain.hpp"
#include "../buffers/prng_buffer.hpp"
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
class Gaugemomentum : public Opencl_Module {
public:
	friend hardware::Device;

	virtual ~Gaugemomentum();

	///////////////////////////////////////////////////
	//Methods on device
	void set_float_to_gaugemomentum_squarenorm_device(const hardware::buffers::Gaugemomentum * in, const hardware::buffers::Plain<hmc_float> * out) const;
	void global_squarenorm_reduction(const hardware::buffers::Plain<hmc_float> * out, const hardware::buffers::Plain<hmc_float> * tmp_buf) const;
	void generate_gaussian_gaugemomenta_device(const hardware::buffers::Gaugemomentum * in, const hardware::buffers::PRNGBuffer * prng) const;
	void set_zero_gaugemomentum(const hardware::buffers::Gaugemomentum *) const;
	/**
	 * This function returns the input gaugemomentum field x
	 * multiplied by alpha and added to the other input gaugemomentum field y
	 * @param x The first input gaugemomentum field (one ae per site)
	 * @param y The second input gaugemomentum field (one ae per site)
	 * @param alpha The real constant
	 * @param out The output gaugemomentum field alpha*x + y (one ae per site)
	 */
	void saxpy_device(const hardware::buffers::Gaugemomentum * x, const hardware::buffers::Gaugemomentum * y, const hardware::buffers::Plain<hmc_float> * alpha, const hardware::buffers::Gaugemomentum * out) const;
	
	/**
	 * Import data from the gaugemomenta array into the given buffer.
	 *
	 * The data in the buffer will be stored in the device specific format.
	 *
	 * @param[out] dest The buffer to write to in the device specific format
	 * @param[in] data The data to write to the buffer
	 */
	void importGaugemomentumBuffer(const hardware::buffers::Gaugemomentum * dest, const ae * const data) const;
	/**
	 * Export data from the given buffer into a normal gaugemomentum array.
	 *
	 * The data in the buffer is assumed to be in the device specific format.
	 *
	 * @param[out] dest An array that the buffer data can be written to.
	 * @param[in] data A buffer containing the data in the device specific format.
	 */
	void exportGaugemomentumBuffer(ae * const dest, const hardware::buffers::Gaugemomentum * buf) const;

	ClSourcePackage get_sources() const noexcept;

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

protected:

	/**
	 * computes work-sizes for a kernel
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
	Gaugemomentum(const hardware::code::OpenClKernelParametersInterface& kernelParams, const hardware::Device * device);

private:
	/**
	 * Collect the kernels for OpenCL.
	 */
	void fill_kernels();
	/**
	 * Clear out the kernels,
	 */
	void clear_kernels();

	ClSourcePackage basic_gaugemomentum_code;

	//kernels
	cl_kernel generate_gaussian_gaugemomenta;
	cl_kernel _set_zero_gaugemomentum;
	cl_kernel gaugemomentum_squarenorm;
	cl_kernel gaugemomentum_squarenorm_reduction;
	cl_kernel gaugemomentum_convert_to_soa;
	cl_kernel gaugemomentum_convert_from_soa;
	cl_kernel gaugemomentum_saxpy;
};

}

}

#endif // _HARDWARE_CODE_GAUGEMOMENTUM_
