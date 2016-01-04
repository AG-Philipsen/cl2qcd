/** @file
 * Complex OpenCL functionality
 *
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.phys.uni-frankfurt.de>
 * Copyright (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
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

#ifndef _HARDWARE_CODE_COMPLEX_
#define _HARDWARE_CODE_COMPLEX_

#include "opencl_module.hpp"

#include "../buffers/plain.hpp"

namespace hardware {

namespace code {

/**
 * An OpenCL device
 *
 * Adds complex operations to basic Opencl_Module class
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Complex : public Opencl_Module {
public:
	friend hardware::Device;

	virtual ~Complex();

	/**
	 * @param in The real number to be transformed in a complex number
	 * @param out A complex number containing "in"
	 */
	void set_complex_to_float_device(const hardware::buffers::Plain<hmc_float> * in, const hardware::buffers::Plain<hmc_complex> * out) const;
	
	/**
	 * @param a The numerator of the fraction
	 * @param b The denominator of the fraction
	 * @param out A complex number containing a/b
	 */
	void set_complex_to_ratio_device(const hardware::buffers::Plain<hmc_complex> * a, const hardware::buffers::Plain<hmc_complex> * b, const hardware::buffers::Plain<hmc_complex> * out) const;
	
	/**
	 * @param a The first term of the multiplication
	 * @param b The second term of the multiplication
	 * @param out A complex number containing a * b
	 */
	void set_complex_to_product_device(const hardware::buffers::Plain<hmc_complex> * a, const hardware::buffers::Plain<hmc_complex> * b, const hardware::buffers::Plain<hmc_complex> * out) const;
	
	/**
	 * @param a The first term of the addition
	 * @param b The second term of the addition
	 * @param out A complex number containing a + b
	 */
	void set_complex_to_sum_device(const hardware::buffers::Plain<hmc_complex> * a, const hardware::buffers::Plain<hmc_complex> * b, const hardware::buffers::Plain<hmc_complex> * out) const;

	/**
	 * @param a The first term of the aubtraction
	 * @param b The second term of the subtraction
	 * @param out A complex number containing a - b
	 */
	void set_complex_to_difference_device(const hardware::buffers::Plain<hmc_complex> * a, const hardware::buffers::Plain<hmc_complex> * b, const hardware::buffers::Plain<hmc_complex> * out) const;
	
	/**
	 * Print the profiling information to a file.
	 *
	 * @param filename Name of file where data is appended.
	 */
	void virtual print_profiling(const std::string& filename, int number) const override;

	ClSourcePackage get_sources() const noexcept;

	/**
	 * Return amount of Floating point operations performed by a specific kernel per call.
	 * NOTE: this is meant to be the "netto" amount in order to be comparable.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual uint64_t get_flop_size(const std::string& in) const override;

	/**
	 * Return amount of bytes read and written by a specific kernel per call.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual size_t get_read_write_size(const std::string& in) const override;

protected:
	/**
	 * Add specific work_size determination for this child class
	 */
	virtual void get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const override;
	/**
	 * @todo: the constructor must be public at the moment in order to be called from OpenClCode class.
	 * 	It may be made private again in the future!
	 */
public:
	Complex(const hardware::code::OpenClKernelParametersInterface& kernelParameters, const hardware::Device * device);

private:
	/**
	 * Collect the kernels for OpenCL.
	 */
	void fill_kernels();

	/**
	 * Clear out the kernels,
	 */
	void clear_kernels();

	ClSourcePackage basic_complex_code;

	//Single
	cl_kernel convert;
	cl_kernel ratio;
	cl_kernel product;
	cl_kernel sum;
	cl_kernel difference;
};

}

}

#endif // _HARDWARE_CODE_COMPLEX_
