/** @file
 * OpenCL device managment for TK kappa and everything executed on them.
 */
#ifndef _MYOPENCL_kH_
#define _MYOPENCL_kH_

#include <cstdlib>
#include <vector>
#include <cstring>
#include <iostream>
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
#include "host_operations_spinor.h"
#include "host_operations_spinorfield.h"
#include "host_operations_fermionmatrix.h"
#include "globaldefs.h"
#include "hmcerrs.h"
#include "types.h"
#include "host_use_timer.h"
#include "host_testing.h"
#include "host_random.h"
#include "opencl.h"

/**
 * An OpenCL device
 *
 * This class wraps all operations on a device. Operations are always specific, e.g. each kernel and copy operation
 * has it's own wrapper function.
 */
class Opencl_k : public opencl {
 
    public:
      virtual hmc_error fill_kernels_file ();
    
//     private:
//         std::vector<std::string> cl_kernels_file;
// 	int isinit;
// 	cl_context context;
// 	cl_command_queue queue;
// 	cl_program clprogram;
// 	cl_kernel heatbath_odd;
// 	cl_kernel heatbath_even;
// 	cl_kernel overrelax_odd;
// 	cl_kernel overrelax_even;
// 	cl_kernel plaquette;
// 	cl_kernel plaquette_reduction;
// 	cl_kernel polyakov;
// 	cl_kernel polyakov_reduction;
// 
// 	//heatbath variables
// 	cl_mem clmem_gaugefield;
// 	cl_mem clmem_rndarray;
// 	cl_mem clmem_plaq;
// 	cl_mem clmem_plaq_buf_glob;
// 	cl_mem clmem_splaq_buf_glob;
// 	cl_mem clmem_tplaq_buf_glob;
// 	cl_mem clmem_splaq;
// 	cl_mem clmem_tplaq;
// 	cl_mem clmem_polyakov;
// 	cl_mem clmem_polyakov_buf_glob;
// 	//!!CP: this is not needed at the moment and since is not copied to the device anywhere!!
// 	cl_mem clmem_theta_gaugefield;
  
};

#endif