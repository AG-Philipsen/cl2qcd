#ifndef _MYOPENCLH_
#define _MYOPENCLH_

#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <CL/cl.hpp>

#include "globaldefs.h"
#include "hmcerrs.h"
#include "types.h"

extern cl_device_type wanted_device;
extern cl_context context;
extern cl_command_queue cmdqueue;
extern cl_program clprogram;
extern std::string cl_kernels_file;

hmc_error init_opencl();
hmc_error finalize_opencl();

#endif
