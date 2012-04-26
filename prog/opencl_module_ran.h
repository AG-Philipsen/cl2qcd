/** @file
 * Heatbath for OpenCL
 */
#ifndef _OPENCLMODULERANH_
#define _OPENCLMODULERANH_

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
#include "inputparameters.h"
#include "opencl_compiler.hpp"

#include "opencl_module.h"

#include "exceptions.h"

/**
 * An OpenCL device
 *
 * Adds random numbers to basic Opencl_Module class
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Opencl_Module_Ran : public Opencl_Module {
public:
	/**
	 * Collect the compiler options for OpenCL.
	 * Virtual method, allows to include more options in inherited classes.
	 */
	virtual void fill_collect_options(stringstream* collect_options);

	/**
	 * Collect the buffers to generate for OpenCL.
	 * Virtual method, allows to include more buffers in inherited classes.
	 */
	virtual void fill_buffers();

	/**
	 * Collect the kernels for OpenCL.
	 * Virtual method, allows to include more kernels in inherited classes.
	 */
	virtual void fill_kernels();

	/**
	 * Clear out the kernels,
	 * Virtual method, allows to clear additional kernels in inherited classes.
	 */
	virtual void clear_kernels();

	/**
	 * Clear out the buffers,
	 * Virtual method, allows to clear additional buffers in inherited classes.
	 */
	virtual void clear_buffers();

	/**
	 * Initialize random arrays.
	 * Called by fill_buffers() so that the arrays are available "in time"
	 */
	void init_random_arrays();

	/**
	 * Get cl_mem object rndarray
	 * @return rndarray
	 */
	cl_mem* get_clmem_rndarray();

	/**
	 * Add specific work_size determination for this child class
	 */
	virtual void get_work_sizes(const cl_kernel kernel, cl_device_type dev_type, size_t * ls, size_t * gs, cl_uint * num_groups);

protected:
	/**
	 * Get number of random states
	 * @return num_rndstates
	 */
	int get_num_rndstates();
	prng_state_dev* rndarray;
	size_t sizeof_rndarray;

private:

	int num_rndstates;
	cl_mem clmem_rndarray;

	/**
	 * Copy the RNG state to the appropriate OpenCL buffer.
	 *
	 * @param host_rndarray The RNG state to copy
	 *         @li HMC_OCLERROR if OpenCL operations fail
	 *         @li HMC_SUCCESS otherwise
	 */
	void copy_rndstate_to_device(prng_state_dev* host_rndarray);

	/**
	 * Copy the RNG state from the OpenCL buffer.
	 *
	 * @param[out] rndarray The RNG copy target
	 *         @li HMC_OCLERROR if OpenCL operations fail
	 *         @li HMC_SUCCESS otherwise
	 */
	void copy_rndstate_from_device(prng_state_dev* rndarray);
};

#endif //OPENCLMODULERANH
