/** @file
 * Heatbath for OpenCL
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
#ifndef _HARDWARE_CODE_KAPPA_
#define _HARDWARE_CODE_KAPPA_

#include "opencl_module.hpp"
#include "../buffers/plain.hpp"
#include "../buffers/su3.hpp"

namespace hardware {

namespace code {

/**
 * An OpenCL device
 *
 * Adds random numbers to basic Opencl_Module class
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Kappa : public Opencl_Module {
public:
	friend hardware::Device;

	virtual ~Kappa();

	/**
	 * Run the calculation of kappa clover. No OpenCL barrier.
	 * @TODO remove beta
	 */
	void run_kappa_clover(const hardware::buffers::Plain<hmc_float> * kappa, const hardware::buffers::SU3 * gaugefield, const hmc_float beta) const;
	/**
	 * Run the calculation of kappa clover. No OpenCL barrier.
	 * @TODO remove beta
	 */
	void run_kappa_clover(const hardware::buffers::SU3 * gaugefield, const hmc_float beta) const;

	/**
	 * Copy kappa_clover from device to host and return it
	 * @return kappa_clover
	 */
	hmc_float get_kappa_clover() const;

protected:
	/**
	 * Return amount of Floating point operations performed by a specific kernel per call.
	 * NOTE: this is meant to be the "netto" amount in order to be comparable.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual uint64_t get_flop_size(const std::string&) const {
		return 0;
	};

	/**
	 * Return amount of bytes read and written by a specific kernel per call.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual size_t get_read_write_size(const std::string&) const {
		return 0;
	};


private:
	/**
	 * Constructor, only to be used by hardware::device
	 *
	 * @param[in] params points to an instance of inputparameters
	 */
	Kappa(const meta::Inputparameters& params, hardware::Device * device);

	/**
	 * Collect the kernels for OpenCL.
	 */
	void fill_kernels();

	/**
	 * Clear out the kernels,
	 */
	void clear_kernels();

	const hardware::buffers::Plain<hmc_float> clmem_kappa_clover;
	cl_kernel kappa_clover_gpu;
};

}

}

#endif // _HARDWARE_CODE_KAPPA_
