/** @file
 * Opencl code working on the gaugefield.
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
#ifndef _HARDWARE_CODE_GAUGEFIELD_
#define _HARDWARE_CODE_GAUGEFIELD_

#include "opencl_module.hpp"
#include "../buffers/su3.hpp"
#include "../buffers/plain.hpp"

namespace hardware {

namespace code {

/**
 * An OpenCL device
 *
 * This class wraps all gaugefield operations on a device. Each kernel
 * has it's own wrapper function.
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Gaugefield : public Opencl_Module {

public:
	friend hardware::Device;

	/**
	 * Destructor.
	 */
	virtual ~Gaugefield();

	// methods which actually calculate something
	/**
	 * Calculate rectangles of a specific gaugefield (on device).
	 *
	 * @deprecated
	 *
	 * @param[in] gf gaugefield to measure on
	 * @param[out] plaq Storage for result of rectangles calculation
	 */
	void gaugeobservables_rectangles(const hardware::buffers::SU3 * gf, hmc_float * const rect) const;
	/**
	 * Calculate plaquette for a specific gaugefield (on device).
	 *
	 * @param[in] gf gaugefield to measure on
	 */
	void plaquette_device(const hardware::buffers::SU3 * gf, const hardware::buffers::Plain<hmc_float> * plaq, const hardware::buffers::Plain<hmc_float> * tplaq, const hardware::buffers::Plain<hmc_float> * splaq) const;
	/**
	 * Calculate rectangles for a specific gaugefield (on device).
	 *
	 * @param[in] gf gaugefield to measure on
	 */
	void rectangles_device(const hardware::buffers::SU3 * gf, const hardware::buffers::Plain<hmc_float> *) const;
	/**
	 * Calculate Polyakov loop for a specific gaugefield (on device).
	 *
	 * @param[in] gf gaugefield to measure on
	 */
	void polyakov_device(const hardware::buffers::SU3 * gf, const hardware::buffers::Plain<hmc_complex> *) const;

	void polyakov_md_local_device(const hardware::buffers::Plain<Matrixsu3> * partial_results, const hardware::buffers::SU3* gf) const;

	void polyakov_md_merge_device(const hardware::buffers::Plain<Matrixsu3> * partial_results, const cl_uint num_slices, const hardware::buffers::Plain<hmc_complex> * pol) const;

	/**
	 * Print the profiling information to a file.
	 *
	 * @param filename Name of file where data is appended.
	 * @param number task-id
	 */
	void virtual print_profiling(const std::string& filename, int number) const override;

	/**
	 * This applies stout smearing to a gaugefield
	 */
	void stout_smear_device(const hardware::buffers::SU3 * in, const hardware::buffers::SU3 * out) const;

	/**
	 * Import the gaugefield data into the OpenCL buffer using the device
	 * specific storage format.
	 *
	 * @param[out] gaugefield The OpenCL buffer to writ the gaugefield data to in the device specific format
	 * @param[in]  data       The gaugefield data to import into the OpenCL buffer
	 *
	 * @todo should not be public
	 */
	void importGaugefield(const hardware::buffers::SU3 * gaugefield, const Matrixsu3 * const data) const;

	/**
	 * Export the gaugefield from the OpenCL buffer, that uses a device
	 * specific storage format, into the given pointer using the generic
	 * storage format.
	 *
	 * @param[out] dest The array to store the gaugefield in
	 */
	void exportGaugefield(Matrixsu3 * const dest, const hardware::buffers::SU3 * gaugefield) const;

	/**
	 * Get the code required to use the gaugefield from kernels.
	 */
	ClSourcePackage get_sources() const noexcept;

protected:
	/**
	 * comutes work-sizes for a kernel
	 * @todo autotune
	 * @param ls local-work-size
	 * @param gs global-work-size
	 * @param num_groups number of work groups
	 * @param dev_type type of device on which the kernel should be executed
	 * @param name name of the kernel for possible autotune-usage, not yet used!!
	 */
	virtual void get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const override;

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

public:
	/**
	 * @param[in] params points to an instance of inputparameters
	 */
	Gaugefield(const hardware::code::OpenClKernelParametersInterface& kernelParams , const hardware::Device * device);
private:

	/**
	 * A set of source files used by all kernels.
	 */
	ClSourcePackage basic_opencl_code;

	//since this is only applicated to the gaugefield, this should be here...
	cl_kernel stout_smear;

	cl_kernel plaquette;
	cl_kernel plaquette_reduction;
	cl_kernel rectangles;
	cl_kernel rectangles_reduction;
	cl_kernel polyakov;
	cl_kernel polyakov_md_local;
	cl_kernel polyakov_md_merge;
	cl_kernel polyakov_reduction;
	cl_kernel convertGaugefieldToSOA;
	cl_kernel convertGaugefieldFromSOA;

	void convertGaugefieldToSOA_device(const hardware::buffers::SU3 * out, const hardware::buffers::Plain<Matrixsu3> * in) const;
	void convertGaugefieldFromSOA_device(const hardware::buffers::Plain<Matrixsu3> * out, const hardware::buffers::SU3 * in) const;

	/**
	 * Collect the kernels for OpenCL.
	 */
	void fill_kernels();
	/**
	 * Clear out the kernels,
	 */
	void clear_kernels();
};

}

}

#endif /* _HARDWARE_CODE_GAUGEFIELD_ */
