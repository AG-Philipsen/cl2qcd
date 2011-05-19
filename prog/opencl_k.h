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
      
	hmc_error run_kappa_karsch_gpu(const hmc_float beta, usetimer * const timer, hmc_float * kappa_karsch_out);
	hmc_error run_kappa_clover_gpu(const hmc_float beta, usetimer * const timer, hmc_float * kappa_clover_out);
	
	cl_mem clmem_kappa_karsch;
	cl_mem clmem_kappa_karsch_buf_glob;
	cl_mem clmem_kappa_clover;
	cl_mem clmem_kappa_clover_buf_glob;
	
	cl_kernel kappa_karsch_gpu;
	cl_kernel kappa_clover_gpu;
 
};

#endif