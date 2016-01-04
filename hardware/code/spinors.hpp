/*
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

/** @file
 * Heatbath for OpenCL
 */
#ifndef _HARDWARE_CODE_SPINORS_
#define _HARDWARE_CODE_SPINORS_

#include "opencl_module.hpp"

#include "../buffers/plain.hpp"
#include "../buffers/spinor.hpp"
#include "../buffers/prng_buffer.hpp"
#include "../size_4.hpp"

namespace hardware {

namespace code {

size_t get_spinorfieldsize(const size_4& params);
size_t get_eoprec_spinorfieldsize(const size_4& params);

/**
 * An OpenCL device
 *
 * Adds random numbers to basic Opencl_Module class
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Spinors : public Opencl_Module {
public:
	friend hardware::Device;
	friend hardware::buffers::Spinor;

	virtual ~Spinors();

	/////////////////////////////////////////
	// device operations
	void generate_gaussian_spinorfield_device(const hardware::buffers::Plain<spinor> * in, const hardware::buffers::PRNGBuffer * prng) const;
	void generate_gaussian_spinorfield_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::PRNGBuffer * prng) const;

	//    linear Algebra operations
	void convert_from_eoprec_device(const hardware::buffers::Spinor * even, const hardware::buffers::Spinor * odd, const hardware::buffers::Plain<spinor> * out) const;
	void convert_to_eoprec_device(const hardware::buffers::Spinor * even, const hardware::buffers::Spinor * odd, const hardware::buffers::Plain<spinor> * in) const;

	void set_complex_to_scalar_product_device(const hardware::buffers::Plain<spinor> * a, const hardware::buffers::Plain<spinor> * b, const hardware::buffers::Plain<hmc_complex> * out) const;
	void set_complex_to_scalar_product_eoprec_device(const hardware::buffers::Spinor * a, const hardware::buffers::Spinor * b, const hardware::buffers::Plain<hmc_complex> * out) const;
	void global_squarenorm_reduction(const hardware::buffers::Plain<hmc_float> * out, const hardware::buffers::Plain<hmc_float> * tmp_buf) const;
	void set_float_to_global_squarenorm_device(const hardware::buffers::Plain<spinor> * a, const hardware::buffers::Plain<hmc_float> * out) const;
	void set_float_to_global_squarenorm_eoprec_device(const hardware::buffers::Spinor * a, const hardware::buffers::Plain<hmc_float> * out) const;
	void set_zero_spinorfield_device(const hardware::buffers::Plain<spinor> * x) const;
	void set_zero_spinorfield_eoprec_device(const hardware::buffers::Spinor * x) const;
	void saxpy_device(const hardware::buffers::Plain<spinor> * x, const hardware::buffers::Plain<spinor> * y, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Plain<spinor> * out) const;
	void saxpy_device(const hardware::buffers::Plain<spinor> * x, const hardware::buffers::Plain<spinor> * y, const hmc_complex alpha, const hardware::buffers::Plain<spinor> * out) const;
	void sax_device(const hardware::buffers::Plain<spinor> * x, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Plain<spinor> * out) const;
	void saxsbypz_device(const hardware::buffers::Plain<spinor> * x, const hardware::buffers::Plain<spinor> * y, const hardware::buffers::Plain<spinor> * z, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Plain<hmc_complex> * beta, const hardware::buffers::Plain<spinor> * out) const;
	void saxpy_eoprec_device(const hardware::buffers::Spinor * x, const hardware::buffers::Spinor * y, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Spinor * out) const;
	void saxpy_eoprec_device(const hardware::buffers::Spinor * x, const hardware::buffers::Spinor * y, const hmc_complex alpha, const hardware::buffers::Spinor * out) const;
	void sax_eoprec_device(const hardware::buffers::Spinor * x, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Spinor * out) const;
	void saxsbypz_eoprec_device(const hardware::buffers::Spinor * x, const hardware::buffers::Spinor * y, const hardware::buffers::Spinor * z, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Plain<hmc_complex> * beta, const hardware::buffers::Spinor * out) const;
	void set_spinorfield_cold_device(const hardware::buffers::Plain<spinor> * inout) const;
	void set_eoprec_spinorfield_cold_device(const hardware::buffers::Spinor * inout) const;

	//    merged kernel calls
	void saxpy_AND_squarenorm_eo_device(const hardware::buffers::Spinor * x, const hardware::buffers::Spinor * y, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Spinor * out, const hardware::buffers::Plain<hmc_complex> * sq_out) const;

	/**
	 * Copy an even-odd preconditioned spinorfield to the given buffer.
	 *
	 * @param buf A buffer of at least get_eoprec_spinorfield_buffer_size() bytes which will
	 *            be filled with the spinorfield in an implementation chosen format.
	 * @param source An array of spinors representing an even-odd field.
	 */
	void copy_to_eoprec_spinorfield_buffer(const hardware::buffers::Spinor * buf, const spinor * const source) const;

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
	Spinors(const hardware::code::OpenClKernelParametersInterface& kernelParameters, const hardware::Device * device);

private:
	/**
	 * Collect the kernels for OpenCL.
	 */
	void fill_kernels();

	/**
	 * Clear out the kernels,
	 */
	void clear_kernels();

	cl_kernel convertSpinorfieldToSOA_eo;
	cl_kernel convertSpinorfieldFromSOA_eo;

	void convertSpinorfieldToSOA_eo_device(const hardware::buffers::Spinor * out, const hardware::buffers::Plain<spinor> * in) const;
	void convertSpinorfieldFromSOA_eo_device(const hardware::buffers::Plain<spinor> * out, const hardware::buffers::Spinor * in) const;

	ClSourcePackage basic_fermion_code;

	//BLAS
	cl_kernel set_spinorfield_cold;
	cl_kernel saxpy;
	cl_kernel saxpy_arg;
	cl_kernel sax;
	cl_kernel saxsbypz;
	cl_kernel set_zero_spinorfield;

	cl_kernel convert_from_eoprec;
	cl_kernel convert_to_eoprec;
	cl_kernel set_eoprec_spinorfield_cold;
	cl_kernel set_zero_spinorfield_eoprec;
	cl_kernel saxpy_eoprec;
	cl_kernel saxpy_arg_eoprec;
	cl_kernel sax_eoprec;
	cl_kernel saxsbypz_eoprec;

	//Scalar Product
	cl_kernel scalar_product;
	cl_kernel scalar_product_reduction;
	cl_kernel global_squarenorm;
	cl_kernel _global_squarenorm_reduction;
	cl_kernel scalar_product_eoprec;
	cl_kernel global_squarenorm_eoprec;

	cl_kernel generate_gaussian_spinorfield;
	cl_kernel generate_gaussian_spinorfield_eo;

	//merged kernels
	cl_kernel saxpy_AND_squarenorm_eo;

	/**
	 * @todo usage of this buffer is dangerous and should probably be semaphored
	 */
	const hardware::buffers::Plain<hmc_complex> * scalar_product_buf;
};

}

}

#endif // _HARDWARE_CODE_SPINORS_
