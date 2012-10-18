/** @file
 * Basic OpenCL functionality
 */
#ifndef _OPENCLMODULEH_
#define _OPENCLMODULEH_

#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>

#include "host_geometry.h"
#include "host_random.h"
#include "host_operations_gaugefield.h"
#include "globaldefs.h"
#include "types.h"
#include "host_use_timer.h"
#include "meta/inputparameters.hpp"
#include "hardware/device.hpp"
#include "opencl_compiler.hpp"
#include "hardware/buffers/su3.hpp"
#include "hardware/buffers/plain.hpp"
#include "meta/util.hpp"

#include "exceptions.h"

/**
 * An OpenCL device
 *
 * This class wraps all operations on a device. Operations are always specific, e.g. each kernel and copy operation
 * has it's own wrapper function.
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Opencl_Module {

public:
	/**
	 * Empty constructor.
	 *
	 * @param[in] params points to an instance of inputparameters
	 */
	Opencl_Module(const meta::Inputparameters& params, hardware::Device * device)
		: parameters(params), device(device), gaugefield(NDIM * meta::get_vol4d(params), device),
		  clmem_rect(1, device), clmem_polyakov(1, device),
		  gf_unsmeared(gaugefield.get_elements(), device),
		  stout_smear(0), rectangles(0), rectangles_reduction(0) {};
	/**
	 * Destructor, calls finalize().
	 *
	 */
	virtual ~Opencl_Module() {
	}

	/**
	 * Free variables. Called by destructor.
	 */
	void finalize();

	/**
	 * Initialize everything. First method to be called.
	 *
	 * @deprecated To be replaced by a proper constructor
	 */
	void init();

	// set and get methods
	/**
	 * Get the context
	 * @return cl_context
	 */
	cl_context get_context();
	/**
	 * Get a pointer to the gaugefield buffer
	 * @return ocl_gaugefield OpenCL buffer with gaugefield
	 */
	const hardware::buffers::SU3 * get_gaugefield();
	/**
	 * Get a pointer to inputparameters
	 * @return parameters
	 */
	const meta::Inputparameters& get_parameters() const noexcept;
	/**
	 * Get OpenCL device
	 * @return device
	 */
	hardware::Device * get_device() const noexcept;

	// methods which actually calculate something
	/**
	 * Calculate plaquette and polyakov of a specific gaugefield (on device).
	 *
	 * @param[out] plaq Storage for result of plaquette calculation
	 * @param[out] tplaq Storage for result of plaquette calculation
	 * @param[out] splaq Storage for result of plaquette calculation
	 * @param[out] pol Storage for result of polyakov calculation
	 */
	void gaugeobservables(hmc_float * const plaq, hmc_float * const tplaq, hmc_float * const splaq, hmc_complex * const pol);
	/**
	 * Calculate plaquette and polyakov of a specific gaugefield (on device).
	 *
	 * @param[in]  gf    The gaugefield on which to compute the observables
	 * @param[out] plaq  Storage for result of plaquette calculation
	 * @param[out] tplaq Storage for result of plaquette calculation
	 * @param[out] splaq Storage for result of plaquette calculation
	 * @param[out] pol   Storage for result of polyakov calculation
	 *
	 * @todo Should not be public
	 */
	void gaugeobservables(const hardware::buffers::SU3 * gf, hmc_float * const plaq, hmc_float * const tplaq, hmc_float * const splaq, hmc_complex * const pol);
	/**
	 * Calculate rectangles of a specific gaugefield (on device).
	 *
	 * @param[in] gf gaugefield to measure on
	 * @param[out] plaq Storage for result of rectangles calculation
	 */
	void gaugeobservables_rectangles(const hardware::buffers::SU3 * gf, hmc_float * const rect);
	/**
	 * Calculate plaquette for a specific gaugefield (on device).
	 *
	 * @param[in] gf gaugefield to measure on
	 */
	void plaquette_device(const hardware::buffers::SU3 * gf, const hardware::buffers::Plain<hmc_float> * plaq, const hardware::buffers::Plain<hmc_float> * tplaq, const hardware::buffers::Plain<hmc_float> * splaq);
	/**
	 * Calculate rectangles for a specific gaugefield (on device).
	 *
	 * @param[in] gf gaugefield to measure on
	 */
	void rectangles_device(const hardware::buffers::SU3 * gf);
	/**
	 * Calculate Polyakov loop for a specific gaugefield (on device).
	 *
	 * @param[in] gf gaugefield to measure on
	 */
	void polyakov_device(const hardware::buffers::SU3 * gf);

	// OpenCL specific methods needed for building/compiling the OpenCL program
	/**
	 * Collect the compiler options for OpenCL.
	 * Virtual method, allows to include more options in inherited classes.
	 */
	virtual void fill_collect_options(std::stringstream* collect_options);
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
	/**
	 * Contains the list of kernel files after call to fill_kernels_file().
	 */
	std::vector<std::string> cl_kernels_file;

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
	 *    clCreateBuffer(context, CL_MEM_ALLOC_HOST_PTR, size, 0, &clerr);
	 *  and returns a pointer to a cl_mem-object located on the host and allocats memory on the host
	 *  if clerr is HMC_SUCCESS
	 *  @param size size of buffer
	 */
	cl_mem create_ahp_buffer(size_t size);

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
	virtual void get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const;

	/**
	 * Return the kernel name as a string
	 * @param[in] kernel
	 * @return kernel_name
	 */
	std::string get_kernel_name(const cl_kernel kernel) const;

	/**
	 * Print the profiling information to a file.
	 *
	 * @param filename Name of file where data is appended.
	 * @param number task-id
	 */
	void virtual print_profiling(const std::string& filename, int number);

	/**
	 * Return amount of Floating point operations performed by a specific kernel per call.
	 * NOTE: this is meant to be the "netto" amount in order to be comparable.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual uint64_t get_flop_size(const std::string& in) const;

	/**
	 * Return amount of bytes read and written by a specific kernel per call.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual size_t get_read_write_size(const std::string& in) const;

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

	/**
	 * This applies stout smearing to a gaugefield
	 */
	void smear_gaugefield(const hardware::buffers::SU3 * gf, const std::vector<const hardware::buffers::SU3 *>& gf_intermediate);
	void stout_smear_device(const hardware::buffers::SU3 * in, const hardware::buffers::SU3 * out);

	/**
	 * This replaces the stout smeared gaugefield with the unsmeared one
	 */
	void unsmear_gaugefield(const hardware::buffers::SU3 * gf);

	usetimer * get_copy_to();
	usetimer * get_copy_on();

	/**
	 * Prints the copy-to/from/on-device-times
	 */
	void print_copy_times(uint64_t totaltime);

	/**
	 * Internal bookeeping function. Only public so it can be called from
	 * C-style callback functions.
	 */
	void markMemReleased(bool host, size_t size);

	/**
	 * Import the gaugefield data into the OpenCL buffer using the device
	 * specific storage format.
	 *
	 * @param[in]  data       The gaugefield data to import into the OpenCL buffer
	 * @param[out] gaugefield The OpenCL buffer to writ the gaugefield data to in the device specific format
	 */
	void importGaugefield(const Matrixsu3 * const data);

	/**
	 * Import the gaugefield data into the OpenCL buffer using the device
	 * specific storage format.
	 *
	 * @param[out] gaugefield The OpenCL buffer to writ the gaugefield data to in the device specific format
	 * @param[in]  data       The gaugefield data to import into the OpenCL buffer
	 *
	 * @todo should not be public
	 */
	void importGaugefield(const hardware::buffers::SU3 * gaugefield, const Matrixsu3 * const data);

	/**
	 * Export the gaugefield from the OpenCL buffer, that uses a device
	 * specific storage format, into the given pointer using the generic
	 * storage format.
	 *
	 * @param[out] dest The array to store the gaugefield in
	 */
	void exportGaugefield(Matrixsu3 * const dest);

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
	TmpClKernel createKernel(const char * const kernel_name, const char * const build_opts = 0);

	const meta::Inputparameters& parameters;

	/**
	 * Print the profiling information for the given kernel to the given file.
	 *
	 * \param filename name of the file to print to
	 * \param kernel the kernel whose information to pring
	 */
	void print_profiling(const std::string& filename, const cl_kernel& kernel) const;

public:

	/**
	 * Get the queue used by this module
	 *
	 * @deprecated queueing shoudl be done via device or buffes
	 */
	cl_command_queue get_queue() const noexcept;

private:

	/**
	 * The device used by this module
	 */
	hardware::Device * const device;

	cl_context ocl_context;

	const hardware::buffers::SU3 gaugefield;

	const hardware::buffers::Plain<hmc_float> clmem_rect;
	const hardware::buffers::Plain<hmc_complex> clmem_polyakov;

	//this is used to save the unsmeared gaugefield if smearing is used
	const hardware::buffers::SU3 gf_unsmeared;

	//since this is only applicated to the gaugefield, this should be here...
	cl_kernel stout_smear;

	cl_kernel plaquette;
	cl_kernel plaquette_reduction;
	cl_kernel rectangles;
	cl_kernel rectangles_reduction;
	cl_kernel polyakov;
	cl_kernel polyakov_reduction;
	cl_kernel convertGaugefieldToSOA;
	cl_kernel convertGaugefieldFromSOA;

	//bunch of timers
	//this is used to measure data-transfer to and from the device
	usetimer copy_to;
	//this is used to measure data-transfer on the device
	usetimer copy_on;

	// memory usage tracing
	size_t allocated_bytes;
	size_t max_allocated_bytes;
	size_t allocated_hostptr_bytes;

	/**
	 * Create an OpenCL buffer object with the given flags,
	 * track memory usage and check errors.
	 */
	cl_mem createBuffer(cl_mem_flags flags, size_t size);

	/**
	 * Create an OpenCL buffer backed by host memory
	 * with the given flags,
	 * track memory usage and check errors.
	 */
	cl_mem createBuffer(cl_mem_flags flags, size_t size, void * host_ptr);

	void convertGaugefieldToSOA_device(const hardware::buffers::SU3 * out, const hardware::buffers::Plain<Matrixsu3> * in);
	void convertGaugefieldFromSOA_device(const hardware::buffers::Plain<Matrixsu3> * out, const hardware::buffers::SU3 * in);
};

#endif //OPENCLMODULEH
