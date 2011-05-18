/** @file
 * OpenCL device managment for TK kappa and everything executed on them.
 */
#ifndef _MYOPENCLKH_
#define _MYOPENCLKH_

#include "opencl.h"

/**
 * An OpenCL device
 *
 * This class wraps all operations on a device. Operations are always specific, e.g. each kernel and copy operation
 * has it's own wrapper function.
 */
class Opencl_k : public Opencl {
 
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