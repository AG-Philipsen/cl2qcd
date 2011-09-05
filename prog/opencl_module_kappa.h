/** @file
 * Heatbath for OpenCL
 */
#ifndef _OPENCLMODULEKAPPAH_
#define _OPENCLMODULEKAPPAH_

#include <cstdlib>
#include <vector>
#include <cstring>
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
class Opencl_Module_Kappa : public Opencl_Module {
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

	void run_kappa_clover(const hmc_float beta);

	cl_mem clmem_kappa_clover;
	cl_mem clmem_kappa_clover_buf_glob;
	cl_kernel kappa_clover_gpu;

protected:

 private:

};

#endif //OPENCLMODULEKAPPAH
