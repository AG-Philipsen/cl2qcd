/** @file
 * Declaration of the hardware::Device class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#ifndef _HARDWARE_DEVICE_HPP_
#define _HARDWARE_DEVICE_HPP_

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include "../meta/inputparameters.hpp"
#include "../opencl_compiler.hpp"

namespace hardware {

	class OptimizationError {
		// thrown to indicate that optimization failed
	};

	/**
	 * Representation of an OpenCL device we are working on.
	 *
	 * Allows the querying and manipulation of the system and
	 * its hardware.
	 */
	class Device {

	public:
		/**
		 * Initialize the device.
		 *
		 * \param context context for the device
		 * \param device_id id of the device to initialize
		 * \param inputparams the input parameters of the application
		 */
		Device(cl_context, cl_device_id, const meta::Inputparameters&);

		~Device();

		/**
		 * Checks whether the device supports double precision.
		 */
		bool is_double_supported();

		// non-copyable
		Device& operator=(const Device&) = delete;
		Device(const Device&) = delete;
		Device() = delete;

		/**
		 * Get the prefered local work size of this device
		 */
		size_t get_preferred_local_thread_num() const noexcept;

		/**
		 * Get the default number of threads on this device
		 */
		size_t get_preferred_global_thread_num() const noexcept;

		/**
		 * Get the number of compute units of the device
		 */
		size_t get_num_compute_units() const noexcept;

		/**
		 * Get the type of the OpenCL device
		 */
		cl_device_type get_device_type() const noexcept;

		/**
		 * Whether this device prefers blocked loops
		 */
		bool get_prefers_blocked_loops() const noexcept;

		/**
		 * Get the name of this device
		 */
		std::string get_name() const noexcept;

		/**
		 * Get the id of this device
		 */
		cl_device_id get_id() const noexcept;

		/**
		 * Create a kernel from source files.
		 *
		 * Usage:
		 * @code
		 * cl_kernel dummy = create_kernel("dummy") << "dummy.cl";
		 * @endcode
		 *
		 * @param kernel_name The name of the kernel to create.
		 * @param build_opts Build options for the kernel
		 */
		TmpClKernel create_kernel(const char * const kernel_name, std::string build_opts) const;

		/**
		 * Enqueue a kernel on the device using the default number of global threads
		 */
		void enqueue_kernel(const cl_kernel kernel) const;

		/**
		 * Enqueue a kernel on the device using the default number of local threads
		 */
		void enqueue_kernel(const cl_kernel kernel, size_t global_threads) const;

		/**
		 * Enqueue a kernel on the device using the default given threads specifications
		 */
		void enqueue_kernel(const cl_kernel kernel, size_t global_threads, size_t local_threads) const;

		/**
		 * Recommend a stride for the given number of elements of the given type
		 *
		 * \param elems      The number of elements to be stored
		 * \param type_size  The size of the basic storage element
		 * \param lane_count How many of the basic storage elements make up the complete type
		 * \return The recommended stride in elements
		 */
		size_t recommend_stride(size_t elems, size_t type_size, size_t lane_count) const;

	private:
		/**
		 * The OpenCL context to be used by this device.
		 */
		const cl_context context;

		/**
		 * The OpenCL device id of the device
		 */
		const cl_device_id device_id;

		/**
		 * The input parameters of the application.
		 */
		const meta::Inputparameters& params;

		/**
		 * The command queue used to perform operations on this device
		 */
		cl_command_queue command_queue;

		/**
		 * The prefered local work size of this device
		 */
		const size_t preferred_local_thread_num;

		/**
		 * The preferred global work size of this device
		 */
		const size_t preferred_global_thread_num;

		/**
		 * The number of compute units of the device
		 */
		const size_t num_compute_units;

		/**
		 * The type of the device
		 */
		const cl_device_type device_type;

		/**
		 * Whether this device supports double precsion
		 */
		const bool supports_double;

		/**
		 * Whether this device prefers blocked loops
		 */
		const bool prefers_blocked_loops;

		/**
		 * The name of this device
		 */
		const std::string name;

		/**
		 * Allow easy use of the command queue
		 */
		operator cl_command_queue() const noexcept;
	};
}

#endif /* _HARDWARE_DEVICE_HPP_ */
