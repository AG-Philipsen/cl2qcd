/** @file
 * OpenCL device managment and everything executed on them.
 */
#ifndef _MYOPENCLH_
#define _MYOPENCLH_

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
#include "hmcerrs.h"
#include "types.h"
#include "host_use_timer.h"
#include "host_random.h"
#include "inputparameters.h"
#include "opencl_compiler.hpp"

/**
 * An OpenCL device
 *
 * This class wraps all operations on a device. Operations are always specific, e.g. each kernel and copy operation
 * has it's own wrapper function.
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Opencl {
public:
	/**
	 * Default constructor that also initializes the device.
	 *
	 * @param wanted The OpenCL device type to be used, e.g. CL_DEVICE_TYPE_CPU or CL_DEVICE_TYPE_GPU
	 * @param timer The timer to use for reporting execution time
	 * @param parameters The parsed input parameters
	 *
	 * @todo Should probably throw an exception on error
	 */
	Opencl(cl_device_type wanted, usetimer* timer, inputparameters* params) {
		this->init(wanted, timer, params);
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
	 * @param timer The timer to use for reporting execution time
	 * @param parameters The parsed input parameters
	 * @return Error code as defined in hmcerrs.h:
	 *         @li HMC_OCLERROR if OpenCL initialization / operations fail
	 *         @li HMC_FILEERROR if one of the kernel files cannot be opened
	 *         @li HMC_SUCCESS otherwise
	 */
	virtual hmc_error init(cl_device_type wanted_device_type, usetimer* timer, inputparameters* parameters);

	/////////////////////////////7
	// communication
	/**
	 * Copy the given gaugefield to the appropriate OpenCL buffer.
	 *
	 * @param host_gaugefield The gaugefield to copy
	 * @param timer The timer to use to measure the copying time
	 * @return Error code as defined in hmcerrs.h:
	 *         @li HMC_OCLERROR if OpenCL operations fail
	 *         @li HMC_SUCCESS otherwise
	 */
	hmc_error copy_gaugefield_to_device(s_gaugefield* gaugefield, usetimer* timer);
//  hmc_error copy_gaugefield_to_device(hmc_gaugefield* host_gaugefield,  usetimer* timer);

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

	/**
	 * Copy the RNG state from the OpenCL buffer.
	 *
	 * @param[out] rndarray The RNG copy target
	 * @param[in,out] timer The timer to use to measure the copying time
	 * @return Error code as defined in hmcerrs.h:
	 *         @li HMC_OCLERROR if OpenCL operations fail
	 *         @li HMC_SUCCESS otherwise
	 */
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
	hmc_error get_gaugefield_from_device(s_gaugefield* gaugefield, usetimer* timer);
//  hmc_error get_gaugefield_from_device(hmc_gaugefield* host_gaugefield,  usetimer* timer);


	////////////////////////////
	//Calculations, calls to kernels
	/**
	 * Perform one heatbath step.
	 */
	hmc_error run_heatbath(const hmc_float beta, usetimer * const timer);

	/**
	 * Perform one overrelaxation step.
	 */
	hmc_error run_overrelax(const hmc_float beta, usetimer * const timer);

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

	/**
	 * Calculate plaquette and polyakov, but with a specific gaugefield.
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
	hmc_error gaugeobservables(cl_mem gf, hmc_float * const plaq, hmc_float * const tplaq, hmc_float * const splaq, hmc_complex * const pol, usetimer * const timer1, usetimer * const timer2);
	/**
	 * returns init status
	 * @return isinit (1==true, 0==false)
	 */
	hmc_error get_init_status();
	/**
	 * Returns private member * parameters
	 * @todo parameters is only used in inherited classes
	 * @return parameters
	 */
	inputparameters * get_parameters ();
	/**
	 * Sets private member * parameters
	 * @return parameters
	 */
	hmc_error set_parameters (inputparameters * parameters_val);

	//protected:

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
	virtual hmc_error fill_kernels(cl_program clprogram);

	///////////////////////////////////////////////
	//get and set methods
	/**
	 * Sets initstatus to 1 (true)
	 *
	 */
	hmc_error set_init_true();
	/**
	 * Sets initstatus to 0 (false)
	 *
	 */
	hmc_error set_init_false();

	//private:

	/**
	 * Called by the destructor.
	 */
	virtual hmc_error finalize();

	/**
	 * Clear out the kernels,
	 * Virtual method, allows to clear additional kernels in inherited classes.
	 */
	virtual hmc_error clear_kernels();

	/**
	 * Clear out the buffers,
	 * Virtual method, allows to clear additional buffers in inherited classes.
	 */
	virtual hmc_error clear_buffers();

	/**
	 * Contains the list of kernel files after call to fill_kernels_file().
	 */
	std::vector<std::string> cl_kernels_file;

	/**
	 * Instance of input_parameters.
	 */
	inputparameters* parameters;

	/** The number of cores (not PEs) of the device */
	cl_uint max_compute_units;

	cl_command_queue queue;

	///////////////////////////////////////////////////////////
	//LZ what follows should eventually be private
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

	cl_platform_id platform;

	/** ID of the OpenCL device wrapped by this object */
	cl_device_id device;
	cl_device_type device_type;
	int isinit;
	cl_context context;
	cl_kernel heatbath_odd;
	cl_kernel heatbath_even;
	cl_kernel overrelax_odd;
	cl_kernel overrelax_even;
	cl_kernel plaquette;
	cl_kernel plaquette_reduction;
	cl_kernel polyakov;
	cl_kernel polyakov_reduction;

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

	/**
	 * Print resource requirements of a kernel object.
	 *
	 * All information is dumped to the trace.
	 *
	 * @param kernel The kernel of which to query the information.
	 */
	void printResourceRequirements(const cl_kernel kernel);

protected:
	/**
	 * A set of source files used by all kernels.
	 */
	ClSourcePackage basic_opencl_code;

	/**
	 * Create a kernel from source files.
	 *
	 * Usage:
	 * @code
	 * cl_kernel dummy = createKernel("dummy") << "dummy.cl";
	 * @endcode
	 *
	 * @param kernel_name The name of the kernel to create.
	 */
	TmpClKernel createKernel(const char * const kernel_name);

private:
	hmc_error init_basic(cl_device_type wanted_device_type, usetimer* timer, inputparameters* parameters);

};

#endif /* _MYOPENCLH_ */
