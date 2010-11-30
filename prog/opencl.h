#ifndef _MYOPENCLH_
#define _MYOPENCLH_

#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <CL/cl.hpp>

#include "globaldefs.h"
#include "hmcerrs.h"
#include "types.h"

extern cl_device_type wanted_device_type;
extern cl_context context;
extern cl_command_queue cmdqueue;
extern cl_program clprogram;
extern std::vector<std::string> const cl_kernels_file;

class opencl {
 public:
  opencl(cl_device_type wanted){init(wanted);};
  ~opencl(){finalize_opencl();};
  hmc_error init(cl_device_type wanted_device_type);
 private:
  //  cl_device_type wanted_device_type;
  cl_context context;
  cl_command_queue cmdqueue;
  cl_program clprogram;
  hmc_error finalize_opencl();
};

#endif
