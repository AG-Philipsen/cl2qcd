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
#include "types.h"
#include "host_use_timer.h"
#include "host_random.h"
#include "inputparameters.h"
#include "opencl_compiler.hpp"

#include "exceptions.h"

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
	 * @param nstates Number of random states
	 *
	 * @todo Should probably throw an exception on error
	 */
	Opencl(cl_device_type wanted, inputparameters* params, int nstates) {
		this->init(wanted, params, nstates);
	};
	/**
	 * Empty constructor. Needed for gaugefield class.
	 */
	Opencl() {};
	~Opencl() {
		finalize();
	};


	/**
	 * Initialize the random array.
	 */
	void init_rndarray(int nstates);

	/**
	 * Initialize the OpenCL device
	 *
	 * @param wanted The OpenCL device type to be used, e.g. CL_DEVICE_TYPE_CPU or CL_DEVICE_TYPE_GPU
	 * @param timer The timer to use for reporting execution time
	 * @param parameters The parsed input parameters
	 *         @li HMC_OCLERROR if OpenCL initialization / operations fail
	 *         @li HMC_FILEERROR if one of the kernel files cannot be opened
	 *         @li HMC_SUCCESS otherwise
	 */
	virtual void init(cl_device_type wanted_device_type, inputparameters* parameters, int nstates);

	/////////////////////////////
	// communication
	/**
	 * Copy the given gaugefield to the appropriate OpenCL buffer.
	 *
	 * @param host_gaugefield The gaugefield to copy
	 * @param timer The timer to use to measure the copying time
	 *         @li HMC_OCLERROR if OpenCL operations fail
	 *         @li HMC_SUCCESS otherwise
	 */
	void copy_gaugefield_to_device(s_gaugefield* gaugefield);
//  void copy_gaugefield_to_device(hmc_gaugefield* host_gaugefield,  usetimer* timer);

	/**
	 * Copy the gaugefield from the device into the given memory location.
	 *
	 * @param host_gaugefield Storage location for the gaugefield
	 *         @li HMC_OCLERROR if OpenCL operations fail
	 *         @li HMC_SUCCESS otherwise
	 */
	void get_gaugefield_from_device(s_gaugefield* gaugefield);
//  void get_gaugefield_from_device(hmc_gaugefield* host_gaugefield,  usetimer* timer);

	/**
	 * Calculate plaquette and polyakov of a specific gaugefield.
	 *
	 * @param[in] gf gaugefield to measure on
	 * @param[out] plaq Storage for result of plaquette calculation
	 * @param[out] tplaq Storage for result of plaquette calculation
	 * @param[out] splaq Storage for result of plaquette calculation
	 * @param[out] pol Storage for result of polyakov calculation
	 */
	void gaugeobservables(cl_mem gf, hmc_float * const plaq, hmc_float * const tplaq, hmc_float * const splaq, hmc_complex * const pol);
	void plaquette_device(cl_mem gf);
	void polyakov_device(cl_mem gf);
	/**
	 * This applies stout smearing to a gaugefield
	 */
	void stout_smear_device();

	/**
	 * returns init status
	 * @return isinit (1==true, 0==false)
	 */
	int get_init_status();
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
	void set_parameters (inputparameters * parameters_val);

	//protected:

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

	///////////////////////////////////////////////
	//get and set methods
	/**
	 * Sets initstatus to 1 (true)
	 *
	 */
	void set_init_true();
	/**
	 * Sets initstatus to 0 (false)
	 *
	 */
	void set_init_false();

	/**
	 * Returns clmem_gaugefield
	 *
	 */
	cl_mem get_clmem_gaugefield();
	/**
	 * Returns device_type
	 */
	cl_device_type get_device_type();

	usetimer * get_copy_to();
	usetimer * get_copy_on();

	//private:

	/**
	 * Called by the destructor.
	 */
	virtual void finalize();

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

	/**
	 *  This calls
	 *    clCreateBuffer(context, CL_MEM_READ_WRITE, size, 0, &clerr)
	 *  and returns a pointer to a read-write cl_mem-object if clerr is HMC_SUCCESS
	 *  @param size size of buffer
	 */
	cl_mem create_rw_buffer(size_t size);

	/**
	 *  This calls
	 *    clCreateBuffer(context, CL_MEM_READ_ONLY, size, 0, &clerr)
	 *  and returns a pointer to a read-only cl_mem-object if clerr is HMC_SUCCESS
	 *  @param size size of buffer
	 */
	cl_mem create_wo_buffer(size_t size);

	/**
	 *  This calls
	 *    clCreateBuffer(context, CL_MEM_WRITE_ONLY, size, 0, &clerr)
	 *  and returns a pointer to a write-only cl_mem-object if clerr is HMC_SUCCESS
	 *  @param size size of buffer
	 */
	cl_mem create_ro_buffer(size_t size);

	/**
	 *  This calls
	 *    clCreateBuffer(context, CL_MEM_USE_HOST_PTR, size, host_pointer, &clerr);
	 *  and returns a pointer to a cl_mem-object located on the host if clerr is HMC_SUCCESS
	 *  @param size size of buffer
	 *  @param host_pointer pointer to memory on host
	 */
	cl_mem create_uhp_buffer(size_t size, void *host_pointer);

	/**
	 *  This calls
	 *    clCreateBuffer(context, CL_MEM_ALLOC_HOST_PTR, size, host_pointer, &clerr);
	 *  and returns a pointer to a cl_mem-object located on the host and allocats memory on the host
	 *  if clerr is HMC_SUCCESS
	 *  @param size size of buffer
	 *  @param host_pointer pointer to memory on host
	 */
	cl_mem create_ahp_buffer(size_t size, void *host_pointer);

	/**
	 *  This calls
	 *    clCreateBuffer(context, CL_MEM_USE_HOST_PTR, size, host_pointer, &clerr);
	 *  and returns a pointer to a cl_mem-object located on the device and
	 *  then copies host-memory pointed to by host-pointer to the device
	 *  if clerr is HMC_SUCCESS
	 *  @param size size of buffer
	 *  @param host_pointer pointer to memory on host
	 */
	cl_mem create_chp_buffer(size_t size, void *host_pointer);

	/**
	 * comutes work-sizes for a kernel
	 * @todo autotune
	 * @param ls local-work-size
	 * @param gs global-work-size
	 * @param num_groups number of work groups
	 * @param dev_type type of device on which the kernel should be executed
	 * @param name name of the kernel for possible autotune-usage, not yet used!!
	 */
	void get_work_sizes(size_t * ls, size_t * gs, cl_uint * num_groups, cl_device_type dev_type, string name = "dummy");
	void get_work_sizes2(const cl_kernel kernel, cl_device_type dev_type, size_t * ls, size_t * gs, cl_uint * num_groups);
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
	string device_double_extension;
	int isinit;
	cl_context context;
	cl_kernel plaquette;
	cl_kernel plaquette_reduction;
	cl_kernel polyakov;
	cl_kernel polyakov_reduction;

	//since this is only applicated to the gaugefield, this should be here...
	cl_kernel stout_smear;

	//bunch of timers
	//this is used to measure data-transfer to and from the device
	usetimer copy_to;
	//this is used to measure data-transfer on the device
	usetimer copy_on;
#ifdef _PROFILING_
	//CP: if PROFILING is activated, one needs a timer for each kernel
	usetimer timer_plaquette;
	usetimer timer_plaquette_reduction;
	usetimer timer_polyakov;
	usetimer timer_polyakov_reduction;

	usetimer timer_stout_smear;
	/**
	 * Return the timer connected to a specific kernel.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual usetimer* get_timer(char * in);

	/**
	   * Return amount of bytes read and written by a specific kernel per call.
	   *
	   * @param in Name of the kernel under consideration.
	   */
	virtual int get_read_write_size(char * in, inputparameters * parameters);

	/**
	 * Print the profiling information of a specific kernel to a file.
	 *
	 * @param filename Name of file where data is appended.
	 * @param kernelName Name of specific kernel.
	 * @param time_total total execution time
	 * @param calls_total total number of kernel calls
	 * @param read_write_size number of bytes read and written by the kernel
	 */
	void print_profiling(std::string filename, char * kernelName, uint64_t time_total, int calls_total, int read_write_size);

	/**
	 * Print the profiling information to a file.
	 *
	 * @param filename Name of file where data is appended.
	 */
	void virtual print_profiling(std::string filename);

#endif

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

	/**
	 * Copy content of a buffer to another buffer inside a queue using
	 *     clEnqueueCopyBuffer(queue, in, out, 0, 0, size , 0, 0, NULL);
	 * @param in source
	 * @param out destination
	 * @param size size of data (out must be equal or bigger than size)
	 */
	void copy_buffer_on_device(cl_mem in, cl_mem out, size_t size);
	/**
	 * Copy content of a buffer on host to a buffer on device inside a queue using
	 *     clEnqueueWriteBuffer(queue, dest, CL_TRUE, 0, size, source, 0, 0, NULL);
	 * This call is a blocking write.
	 * @param source
	 * @param dest
	 * @param size size of data (out must be equal or bigger than size)
	 */
	void copy_buffer_to_device(void * source, cl_mem dest, size_t size);
	/**
	 * Copy content of a buffer on device to a buffer on host inside a queue using
	 *    clEnqueueReadBuffer(queue, source, CL_TRUE, 0, size, dest, 0, NULL, NULL);
	 * This call is a blocking read.
	 * @param source
	 * @param dest
	 * @param size size of data (out must be equal or bigger than size)
	 */
	void get_buffer_from_device(cl_mem source, void * dest, size_t size);

	int get_num_rndstates();

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

	/**
	 * Get number of threads
	 * @return int numthreads
	 *
	 */
	int get_numthreads();

private:
	void init_basic(cl_device_type wanted_device_type, inputparameters* parameters, int nstates);

	int num_rndstates;

	int numthreads;

};

#endif /* _MYOPENCLH_ */
