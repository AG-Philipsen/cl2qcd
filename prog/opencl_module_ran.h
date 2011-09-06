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
#include "host_operations_complex.h"
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

  void init_random_arrays();


  hmc_ocl_ran* rndarray;
  size_t sizeof_rndarray;

	/**
	 * Copy the RNG state to the appropriate OpenCL buffer.
	 *
	 * @param host_rndarray The RNG state to copy
	 *         @li HMC_OCLERROR if OpenCL operations fail
	 *         @li HMC_SUCCESS otherwise
	 */
	void copy_rndarray_to_device(hmc_ocl_ran* host_rndarray);

	/**
	 * Copy the RNG state from the OpenCL buffer.
	 *
	 * @param[out] rndarray The RNG copy target
	 *         @li HMC_OCLERROR if OpenCL operations fail
	 *         @li HMC_SUCCESS otherwise
	 */
	void copy_rndarray_from_device(hmc_ocl_ran* rndarray);	

	cl_mem* get_clmem_rndarray();

protected:
	int get_num_rndstates();

 private:

	int num_rndstates;
	cl_mem clmem_rndarray;

};

#endif //OPENCLMODULERANH
