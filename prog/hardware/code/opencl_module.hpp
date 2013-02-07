/** @file
 * Functionality common to all Opencl_Modules
 * Basically just a utility for keeping the device and input parameters references as well as some profiling utilities.
 */

#ifndef _HARDWARE_CODE_OPENCLMODULE_
#define _HARDWARE_CODE_OPENCLMODULE_

#include <string>

#include "../../meta/inputparameters.hpp"
#include "../../opencl_compiler.hpp"

// predeclaration as headers only use pointers and friend to this
namespace hardware {
class Device;

namespace code {

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
	 * Get a pointer to inputparameters
	 * @return parameters
	 */
	const meta::Inputparameters& get_parameters() const noexcept;

	/**
	 * Get OpenCL device
	 * @return device
	 */
	hardware::Device * get_device() const noexcept;

	/**
	 * Print the profiling information to a file.
	 *
	 * @param filename Name of file where data is appended.
	 * @param number task-id
	 */
	void virtual print_profiling(const std::string& filename, int number) const;

protected:
	/**
	 * Protected constructor to keep this class abstract.
	 *
	 * @param[in] params points to an instance of inputparameters
	 */
	Opencl_Module(const meta::Inputparameters& params, hardware::Device * device)
		: parameters(params), device(device) { };

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
	 * Return amount of Floating point operations performed by a specific kernel per call.
	 * NOTE: this is meant to be the "netto" amount in order to be comparable.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual uint64_t get_flop_size(const std::string& in) const = 0;

	/**
	 * Return amount of bytes read and written by a specific kernel per call.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual size_t get_read_write_size(const std::string& in) const = 0;

	/**
	 * Return the kernel name as a string
	 * @param[in] kernel
	 * @return kernel_name
	 */
	std::string get_kernel_name(const cl_kernel kernel) const;

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
	TmpClKernel createKernel(const char * const kernel_name, std::string build_opts = "") const;

	/**
	 * Print the profiling information for the given kernel to the given file.
	 *
	 * \param filename name of the file to print to
	 * \param kernel the kernel whose information to pring
	 */
	void print_profiling(const std::string& filename, const cl_kernel& kernel) const;

private:

	/**
	 * The input parameters this modules is parametrized for.
	 */
	const meta::Inputparameters& parameters;

	/**
	 * The device used by this module
	 */
	hardware::Device * const device;
};

}

}

#endif //OPENCLMODULEH
