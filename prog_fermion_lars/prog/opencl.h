#ifndef _MYOPENCLH_
#define _MYOPENCLH_

#include <cstdlib>
#include <vector>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <CL/cl.h>

#include "operations.h"
#include "globaldefs.h"
#include "hmcerrs.h"
#include "types.h"

//give a list of all kernel-files
std::vector<std::string> const cl_kernels_file = {"opencl_header.cl","operations_kernels.cl","opencl_heatbath.cl","testing_kernel.cl"};



class opencl {
 public:
  opencl(cl_device_type wanted){init(wanted);};
  ~opencl(){finalize();};
  hmc_error init(cl_device_type wanted_device_type);
  hmc_error copy_gaugefield_to_device(hmc_gaugefield* host_gaugefield);
  hmc_error get_gaugefield_from_device(hmc_gaugefield* host_gaugefield);
  hmc_error run_heatbath(int nsteps, double beta);
  hmc_error testing();
  hmc_error finalize();
 private:
  int isinit;
  cl_context context;
  cl_command_queue queue;
  cl_program clprogram;
  cl_mem clmem_gaugefield;
  cl_mem clmem_beta;
  cl_kernel heatbath;
};

#endif
