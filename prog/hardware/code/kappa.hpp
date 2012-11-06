/** @file
 * Heatbath for OpenCL
 */
#ifndef _HARDWARE_CODE_KAPPA_
#define _HARDWARE_CODE_KAPPA_

#include "../../types.h"

#include "opencl_module.hpp"
#include "../buffers/plain.hpp"
#include "../buffers/su3.hpp"

/**
 * An OpenCL device
 *
 * Adds random numbers to basic Opencl_Module class
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Opencl_Module_Kappa : public hardware::code::Opencl_Module {
public:
	friend hardware::Device;

	virtual ~Opencl_Module_Kappa();

	/**
	 * Run the calculation of kappa clover. No OpenCL barrier.
	 * @TODO remove beta
	 */
	void run_kappa_clover(const hardware::buffers::SU3 * gaugefield, const hmc_float beta);

	/**
	 * Copy kappa_clover from device to host and return it
	 * @return kappa_clover
	 */
	hmc_float get_kappa_clover();

protected:
	/**
	 * Return amount of Floating point operations performed by a specific kernel per call.
	 * NOTE: this is meant to be the "netto" amount in order to be comparable.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual uint64_t get_flop_size(const std::string&) const { return 0; };

	/**
	 * Return amount of bytes read and written by a specific kernel per call.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual size_t get_read_write_size(const std::string&) const { return 0; };


private:
	/**
	 * Constructor, only to be used by hardware::device
	 *
	 * @param[in] params points to an instance of inputparameters
	 */
	Opencl_Module_Kappa(const meta::Inputparameters& params, hardware::Device * device);

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

#endif // _HARDWARE_CODE_KAPPA_
