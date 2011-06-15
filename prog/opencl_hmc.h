/** @file
 * OpenCL device managment and everything executed on them -- including HMC-algorithm.
 */
#ifndef _MYOPENCLHMCH_
#define _MYOPENCLHMCH_

//these are the headers of the mother-classes
#include "opencl.h"
#include "opencl_fermions.h"
//CP: this includes the struct-definitions for the spinors...
#include "types_hmc.h"

/**
 * An OpenCL device for the HMC algorithm.
 *
 * This class wraps all operations on a device. Inherited from classes Opencl and Opencl_fermions.
 */
class Opencl_hmc : public Opencl_fermions {
  public:
	 /**
	 * Collect a vector of kernel file names.
	 * Virtual method, allows to include more kernel files in inherited classes.
	 */
	virtual hmc_error fill_kernels_file ();
	/**
	 * Collect the compiler options for OpenCL.
	 * Virtual method, allows to include more options in inherited classes.
	 */
	virtual hmc_error fill_collect_options(stringstream* collect_options);
	/**
	 * Collect the buffers to generate for OpenCL.
	 * Virtual method, allows to include more buffers in inherited classes.
	 */
	virtual hmc_error fill_buffers();
	/**
	 * Collect the kernels for OpenCL.
	 * Virtual method, allows to include more kernels in inherited classes.
	 */
	virtual hmc_error fill_kernels();


	/**
	 * Initialize the OpenCL device including fermion capabilities
	 *
	 * @param wanted The OpenCL device type to be used, e.g. CL_DEVICE_TYPE_CPU or CL_DEVICE_TYPE_GPU
	 * @param timer The timer to use for reporting execution time
	 * @param parameters The parsed input parameters
	 * @return Error code as defined in hmcerrs.h:
	 *         @li HMC_OCLERROR if OpenCL initialization / operations fail
	 *         @li HMC_FILEERROR if one of the kernel files cannot be opened
	 *         @li HMC_SUCCESS otherwise
	 */
	virtual hmc_error init(cl_device_type wanted_device_type, usetimer* timer, inputparameters* parameters);

	hmc_error finalize_hmc();

	hmc_error init_hmc_variables(inputparameters* parameters, usetimer* timer);

};
#endif // _MYOPENCLHMCH_
