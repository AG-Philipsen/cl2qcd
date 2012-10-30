/** @file
 * Heatbath for OpenCL
 */
#ifndef _OPENCLMODULEKAPPAH_
#define _OPENCLMODULEKAPPAH_

#include <cstdlib>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include "host_geometry.h"
#include "host_operations_gaugefield.h"
#include "globaldefs.h"
#include "types.h"
#include "host_use_timer.h"
#include "host_random.h"
#include "opencl_compiler.hpp"

#include "opencl_module_gaugefield.h"

#include "exceptions.h"

/**
 * An OpenCL device
 *
 * Adds random numbers to basic Opencl_Module class
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Opencl_Module_Kappa : public Opencl_Module {
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

#endif //OPENCLMODULEKAPPAH
