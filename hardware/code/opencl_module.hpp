/** @file
 * Functionality common to all Opencl_Modules
 * Basically just a utility for keeping the device and input parameters references as well as some profiling utilities.
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _HARDWARE_CODE_OPENCLMODULE_
#define _HARDWARE_CODE_OPENCLMODULE_

#include "../../common_header_files/types.h"

#include "../openClKernelParameters.hpp"
#include "../opencl_compiler.hpp"
#include "../../host_functionality/logger.hpp"
#include "../device.hpp"

#include "../buffers/3x3.hpp"
#include "../buffers/su3.hpp"
#include "../buffers/prng_buffer.hpp"
#include "../buffers/spinor.hpp"
#include "../buffers/su3vec.hpp"
#include "../buffers/gaugemomentum.hpp"

#include <string>
#include <limits>
#include <fstream>
#include <cmath>

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
	 * Get OpenCL device
	 * @return device
	 */
	const hardware::Device * get_device() const noexcept;

	/**
	 * Print the profiling information to a file.
	 *
	 * @param filename Name of file where data is appended.
	 * @param number task-id
	 */
	void virtual print_profiling(const std::string& filename, int number) const;
	
	/**
	 * Returns the sources for all children modules.
	 * @return basic_sources
	 */
	ClSourcePackage get_basic_sources() const noexcept;

protected:
	Opencl_Module(const hardware::code::OpenClKernelParametersInterface& kernelParameters, const hardware::Device * device);
	~Opencl_Module();

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

	const hardware::code::OpenClKernelParametersInterface * kernelParameters;

private:

	/**
	 * The device used by this module
	 */
	const hardware::Device * device;
	
	/**
	 * The basic source used by all modules (children classes of this)
	 */
	ClSourcePackage basic_sources;

};

template<typename T> T module_metric_not_implemented() {
	return std::numeric_limits<T>::max();
}

}

}

#endif //OPENCLMODULEH
