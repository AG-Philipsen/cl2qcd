/** @file
 * OpenCL device managment and everything executed on them.
 */
#ifndef _MYOPENCLH_
#define _MYOPENCLH_

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

/**
 * An OpenCL device
 *
 * This class wraps all operations on a device. Operations are always specific, e.g. each kernel and copy operation
 * has it's own wrapper function.
 */
class Opencl {
public:
	/**
	 * Default constructor that also initializes the device.
	 *
	 * @param wanted The OpenCL device type to be used, e.g. CL_DEVICE_TYPE_CPU or CL_DEVICE_TYPE_GPU
	 * @param ls The local work size to be used on the device
	 * @param gs The global work size to be used on the device
	 * @param timer The timer to use for reporting execution time
	 * @param parameters The parsed input parameters
	 *
	 * @todo Should probably throw an exception on error
	 */
	Opencl(cl_device_type wanted, const size_t ls, const size_t gs, usetimer* timer, inputparameters* parameters) {
		init(wanted, ls, gs, timer, parameters);
	};
	/**
	 * Empty constructor. Needed for gaugefield class.
	 */
	Opencl() {};
	~Opencl() {
		finalize();
	};

	/**
	 * Initialize the OpenCL device
	 *
	 * @param wanted The OpenCL device type to be used, e.g. CL_DEVICE_TYPE_CPU or CL_DEVICE_TYPE_GPU
	 * @param ls The local work size to be used on the device
	 * @param gs The global work size to be used on the device
	 * @param timer The timer to use for reporting execution time
	 * @param parameters The parsed input parameters
	 * @return Error code as defined in hmcerrs.h:
	 *         @li HMC_OCLERROR if OpenCL initialization / operations fail
	 *         @li HMC_FILEERROR if one of the kernel files cannot be opened
	 *         @li HMC_SUCCESS otherwise
	 */
	hmc_error init(cl_device_type wanted_device_type, const size_t local_work_size, const size_t global_work_size, usetimer* timer, inputparameters* parameters);

	/**
	 * Copy the given gaugefield to the appropriate OpenCL buffer.
	 *
	 * @param host_gaugefield The gaugefield to copy
	 * @param timer The timer to use to measure the copying time
	 * @return Error code as defined in hmcerrs.h:
	 *         @li HMC_OCLERROR if OpenCL operations fail
	 *         @li HMC_SUCCESS otherwise
	 */
	hmc_error copy_gaugefield_to_device(hmc_gaugefield* host_gaugefield,  usetimer* timer);

	/**
	 * Copy the RNG state to the appropriate OpenCL buffer.
	 *
	 * @param host_rndarray The RNG state to copy
	 * @param timer The timer to use to measure the copying time
	 * @return Error code as defined in hmcerrs.h:
	 *         @li HMC_OCLERROR if OpenCL operations fail
	 *         @li HMC_SUCCESS otherwise
	 */
	hmc_error copy_rndarray_to_device(hmc_rndarray host_rndarray,  usetimer* timer);

	hmc_error copy_rndarray_from_device(hmc_rndarray rndarray, usetimer* timer);

	/**
	 * Copy the gaugefield from the device into the given memory location.
	 *
	 * @param host_gaugefield Storage location for the gaugefield
	 * @param timer The timer to use to measure the copying time
	 * @return Error code as defined in hmcerrs.h:
	 *         @li HMC_OCLERROR if OpenCL operations fail
	 *         @li HMC_SUCCESS otherwise
	 */
	hmc_error get_gaugefield_from_device(hmc_gaugefield* host_gaugefield,  usetimer* timer);

	/**
	 * Perform one heatbath step.
	 */
	hmc_error run_heatbath(const hmc_float beta, usetimer * const timer);

	/**
	 * Perform one overrelaxation step.
	 */
	hmc_error run_overrelax(const hmc_float beta, usetimer * const timer);

	/**
	 * Compute the transport coefficient kappa with the energy-momentum-tensor discretized by Karsch&Wyld on GPU
	 * @param[in] local_work_size OpenCL local_work_size
	 * @param[in] global_work_size OpenCL global_work_size
	 * @param[in,out] timer time for measurement
	 * @return Error code as defined in hmcerrs.h
	 */	
	hmc_error run_kappa_karsch_gpu(const size_t local_work_size, const size_t global_work_size, usetimer* timer_karsch);
	
	/**
	 * Compute the transport coefficient kappa with the energy-momentum-tensor discretized by Karsch&Wyld  on GPU
	 * @param[in] local_work_size OpenCL local_work_size
	 * @param[in] global_work_size OpenCL global_work_size
	 * @param[in,out] timer time for measurement
	 * @return Error code as defined in hmcerrs.h
	 */	
	hmc_error run_kappa_clover_gpu (const size_t local_work_size, const size_t global_work_size, usetimer* timer_clover);
  
	/**
	 * Calculate plaquette and polyakov.
	 *
	 * @param[out] plaq Storage for result of plaquette calculation
	 * @param[out] tplaq Storage for result of plaquette calculation
	 * @param[out] splaq Storage for result of plaquette calculation
	 * @param[out] pol Storage for result of polyakov calculation
	 * @param[in,out] timer1 Timer into which to aggregate plaquette calculation time
	 * @param[in,out] timer2 Timer into which to aggregate polyakov calculation time
	 * @return Error code as defined in hmcerrs.h:
	 *         @li HMC_OCLERROR if OpenCL operations fail
	 *         @li HMC_SUCCESS otherwise
	 */
	hmc_error gaugeobservables(hmc_float * const plaq, hmc_float * const tplaq, hmc_float * const splaq, hmc_complex * const pol, usetimer * const timer1, usetimer * const timer2);

	virtual hmc_error fill_kernels_file ();
	

	hmc_error finalize();
// private:
        std::vector<std::string> cl_kernels_file;
	int isinit;
	cl_context context;
	cl_command_queue queue;
	cl_program clprogram;
	cl_kernel heatbath_odd;
	cl_kernel heatbath_even;
	cl_kernel overrelax_odd;
	cl_kernel overrelax_even;
	cl_kernel plaquette;
	cl_kernel plaquette_reduction;
	cl_kernel polyakov;
	cl_kernel polyakov_reduction;

	//heatbath variables
	cl_mem clmem_gaugefield;
	cl_mem clmem_rndarray;
	cl_mem clmem_plaq;
	cl_mem clmem_plaq_buf_glob;
	cl_mem clmem_splaq_buf_glob;
	cl_mem clmem_tplaq_buf_glob;
	cl_mem clmem_splaq;
	cl_mem clmem_tplaq;
	cl_mem clmem_polyakov;
	cl_mem clmem_polyakov_buf_glob;
	//!!CP: this is not needed at the moment and since is not copied to the device anywhere!!
	cl_mem clmem_theta_gaugefield;


	/** The number of cores (not PEs) of the device */
	cl_uint max_compute_units;

	/**
	 * Enqueue the given kernel on the device. Local work size will be determined
	 * automatically from device and kernel properties.
	 *
	 * @param kernel The kernel to execute.
	 * @param global_work_size The number of threads to run.
	 *
	 * @todo local work size decision might need ot become less automatic
	 * @todo global work size will also depend on device ...
	 */
	void enqueueKernel(const cl_kernel kernel, const size_t global_work_size);

	/**
	 * Enqueue the given kernel on the device. Local work size will be determined
	 * automatically from device and kernel properties.
	 *
	 * @param kernel The kernel to execute.
	 * @param global_work_size The number of threads to run.
	 *
	 * @todo local work size decision might need ot become less automatic
	 * @todo global work size will also depend on device ...
	 */
	void enqueueKernel(const cl_kernel kernel, const size_t global_work_size, const size_t local_work_size);
};

#endif /* _MYOPENCLH_ */
