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
	void rectangles_device(const hardware::buffers::SU3 * gf, const hardware::buffers::Plain<hmc_float> *);
	/**
	 * Calculate Polyakov loop for a specific gaugefield (on device).
	 *
	 * @param[in] gf gaugefield to measure on
	 */
	void polyakov_device(const hardware::buffers::SU3 * gf, const hardware::buffers::Plain<hmc_complex> *);

	// OpenCL specific methods needed for building/compiling the OpenCL program
	/**
	 * Collect the compiler options for OpenCL.
	 * Virtual method, allows to include more options in inherited classes.
	 */
	virtual void fill_collect_options(std::stringstream* collect_options);
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
	 * Contains the list of kernel files after call to fill_kernels_file().
	 */
	std::vector<std::string> cl_kernels_file;

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
	 * This applies stout smearing to a gaugefield
	 */
	void smear_gaugefield(const hardware::buffers::SU3 * gf, const std::vector<const hardware::buffers::SU3 *>& gf_intermediate);
	void stout_smear_device(const hardware::buffers::SU3 * in, const hardware::buffers::SU3 * out);

	/**
	 * This replaces the stout smeared gaugefield with the unsmeared one
	 */
	void unsmear_gaugefield(const hardware::buffers::SU3 * gf);

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

	const meta::Inputparameters& parameters;

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

	/**
	 * Print the profiling information for the given kernel to the given file.
	 *
	 * \param filename name of the file to print to
	 * \param kernel the kernel whose information to pring
	 */
	void print_profiling(const std::string& filename, const cl_kernel& kernel) const;

private:

	/**
	 * The device used by this module
	 */
	hardware::Device * const device;

	const hardware::buffers::SU3 gaugefield;

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

	void convertGaugefieldToSOA_device(const hardware::buffers::SU3 * out, const hardware::buffers::Plain<Matrixsu3> * in);
	void convertGaugefieldFromSOA_device(const hardware::buffers::Plain<Matrixsu3> * out, const hardware::buffers::SU3 * in);
};

#endif //OPENCLMODULEH
