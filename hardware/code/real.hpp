/** @file
 * Real OpenCL functionality
 *
 * Copyright (c) 2014 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
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

#ifndef _HARDWARE_CODE_REAL_
#define _HARDWARE_CODE_REAL_

#include "opencl_module.hpp"

#include "../buffers/plain.hpp"

namespace hardware {

namespace code {

/**
 * An OpenCL device
 *
 * This module contains algebra operations with real numbers that are needed somewhere in the code
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Real : public Opencl_Module {
public:
	friend hardware::Device;

	virtual ~Real();
	
	/**
	 * Function to get one element of a buffer
	 */
	void set_real_to_vector_element_device(const hardware::buffers::Plain<hmc_float> * in, const int index, const hardware::buffers::Plain<hmc_float> * out) const;
	
	/**
	 * Function to set one element of a buffer
	 */
	void set_vector_element_to_real_device(const hardware::buffers::Plain<hmc_float> * in, const int index, const hardware::buffers::Plain<hmc_float> * out) const;
	
	/**
	 * @param a The numerator of the fraction
	 * @param b The denominator of the fraction
	 * @param out A real number containing a/b
	 */
	void set_real_to_ratio_device(const hardware::buffers::Plain<hmc_float> * a, const hardware::buffers::Plain<hmc_float> * b, const hardware::buffers::Plain<hmc_float> * out) const;
	
	/**
	 * @param a The first term of the multiplication
	 * @param b The second term of the multiplication
	 * @param out A real number containing a * b
	 */
	void set_real_to_product_device(const hardware::buffers::Plain<hmc_float> * a, const hardware::buffers::Plain<hmc_float> * b, const hardware::buffers::Plain<hmc_float> * out) const;
	
	/**
	 * @param a The first term of the addition
	 * @param b The second term of the addition
	 * @param out A real number containing a + b
	 */
	void set_real_to_sum_device(const hardware::buffers::Plain<hmc_float> * a, const hardware::buffers::Plain<hmc_float> * b, const hardware::buffers::Plain<hmc_float> * out) const;

	/**
	 * @param a The first term of the aubtraction
	 * @param b The second term of the subtraction
	 * @param out A real number containing a - b
	 */
	void set_real_to_difference_device(const hardware::buffers::Plain<hmc_float> * a, const hardware::buffers::Plain<hmc_float> * b, const hardware::buffers::Plain<hmc_float> * out) const;

	/**
	 * Tool for the multimass conjugate gradient algorithm.
	 * @param zeta_prev zeta at the previous iteration (vector of real numbers)
	 * @param zeta_prev_prev zeta at the previous but one iteration (vector of real numbers)
	 * @param sbeta_prev scalar beta at the previous iteration (real number)
	 * @param sbeta_pres scalar beta at the present iteration (real number)
	 * @param salpha_prev scalar alpha at the previous iteration (real number)
	 * @param sigma shifts of the cgm (vector of real numbers)
	 * @param numeq number of masses of the cgm (unsigned integer)
	 * @param out (zeta_prev_prev[k] * zeta_prev[k] * sbeta_prev) / 
	 *            / (sbeta_pres * salpha_prev * (zeta_prev_prev[k] - zeta_prev[k]) +
	 *               + zeta_prev_prev[k] * sbeta_prev * (1 - sigma[k] * beta_pres))
	 *            
	 *            vector of real numbers.
	 * @note pres means updated in the iteration in progress
	 *       prev means updated in the previous iteration
	 *       prev_prev means updated in the previous but one iteration
	 */
	void update_zeta_cgm_device(const hardware::buffers::Plain<hmc_float> * zeta_prev, const hardware::buffers::Plain<hmc_float> * zeta_prev_prev, const hardware::buffers::Plain<hmc_float> * sbeta_prev, const hardware::buffers::Plain<hmc_float> * sbeta_pres, const hardware::buffers::Plain<hmc_float> * salpha_prev, const hardware::buffers::Plain<hmc_float> * sigma, const int numeq, const hardware::buffers::Plain<hmc_float> * out) const;
	
	/**
	 * Tool for the multimass conjugate gradient algorithm.
	 * @param sbeta_pres scalar beta at the present iteration (real number)
	 * @param zeta_pres zeta at the present iteration (vector of real numbers) 
	 * @param zeta_prev zeta at the previous iteration (vector of real numbers)
	 * @param numeq number of masses of the cgm (unsigned integer)
	 * @param out sbeta_prev * zeta_pres[k] / zeta_prev[k]   (vector of real numbers)
	 */
	void update_beta_cgm_device(const hardware::buffers::Plain<hmc_float> * sbeta_pres, const hardware::buffers::Plain<hmc_float> * zeta_pres, const hardware::buffers::Plain<hmc_float> * zeta_prev, const int numeq, const hardware::buffers::Plain<hmc_float> * out) const;
	
	/**
	 * Tool for the multimass conjugate gradient algorithm.
	 * @param salpha_pres scalar alpha at the present iteration (real number)
	 * @param zeta_pres zeta at the present iteration (vector of real numbers) 
	 * @param beta_pres beta at the present iteration (vector of real numbers)
	 * @param zeta_prev zeta at the previous iteration (vector of real numbers) 
	 * @param sbeta_pres scalar beta at the present iteration (real number)
	 * @param numeq number of masses of the cgm (unsigned integer)
	 * @param out salpha_pres * zeta_pres[k] * beta_pres[k] / (zeta_prev[k] * sbeta_pres)  (vector of real numbers)
	 */
	void update_alpha_cgm_device(const hardware::buffers::Plain<hmc_float> * salpha_pres, const hardware::buffers::Plain<hmc_float> * zeta_pres, const hardware::buffers::Plain<hmc_float> * beta_pres, const hardware::buffers::Plain<hmc_float> * zeta_prev, const hardware::buffers::Plain<hmc_float> * sbeta_pres, const int numeq, const hardware::buffers::Plain<hmc_float> * out) const;
	
	  
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
	 * @param numeq This is the parameter numeq of the update_***_cgm_device methods.
	 * @todo Think how avoid this second parameter
	 */
	virtual uint64_t get_flop_size_update(const std::string& in, const int numeq) const;
	virtual uint64_t get_flop_size(const std::string& in) const override;

	/**
	 * Return amount of bytes read and written by a specific kernel per call.
	 *
	 * @param in Name of the kernel under consideration.
	 * @param numeq This is the parameter numeq of the update_***_cgm_device methods.
	 * @todo Think how avoid this second parameter
	 */
	virtual size_t get_read_write_size_update(const std::string& in, const int numeq) const;
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
	Real(const hardware::code::OpenClKernelParametersInterface& kernelParameters, const hardware::Device * device);
private:
	/**
	 * Collect the kernels for OpenCL.
	 */
	void fill_kernels();

	/**
	 * Clear out the kernels,
	 */
	void clear_kernels();

	ClSourcePackage basic_real_code;

	//Setting operations
	cl_kernel get_elem_vec;
	cl_kernel set_elem_vec;
	
	//Single operations
	cl_kernel ratio;
	cl_kernel product;
	cl_kernel sum;
	cl_kernel difference;
	
	//Update cgm kernels
	cl_kernel update_alpha_cgm;
	cl_kernel update_beta_cgm;
	cl_kernel update_zeta_cgm;

};

}

}

#endif // _HARDWARE_CODE_REAL_
