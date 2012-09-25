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

#include <map>
#include "../meta/inputparameters.hpp"
#include "../opencl_compiler.hpp"
#include "profiling_data.hpp"

/**
 * Forward dekleration of friends
 */
class Opencl_Module;

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

		friend ::Opencl_Module; // during refactoring this needs access to the ocl-queue

	public:
		/**
		 * Initialize the device.
		 *
		 * \param context context for the device
		 * \param device_id id of the device to initialize
		 * \param inputparams the input parameters of the application
		 * \param enable_profiling enable profiling on this device
		 */
		Device(cl_context, cl_device_id, const meta::Inputparameters&, bool enable_profiling = false);

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
		 * Whether this device prefers SOA storage
		 *
		 * @todo This should be datatype specific
		 */
		bool get_prefers_soa() const noexcept;

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
		void enqueue_kernel(const cl_kernel kernel);

		/**
		 * Enqueue a kernel on the device using the default number of local threads
		 */
		void enqueue_kernel(const cl_kernel kernel, size_t global_threads);

		/**
		 * Enqueue a kernel on the device using the default given threads specifications
		 */
		void enqueue_kernel(const cl_kernel kernel, size_t global_threads, size_t local_threads);

		/**
		 * Recommend a stride for the given number of elements of the given type
		 *
		 * \param elems      The number of elements to be stored
		 * \param type_size  The size of the basic storage element
		 * \param lane_count How many of the basic storage elements make up the complete type
		 * \return The recommended stride in elements
		 */
		size_t recommend_stride(size_t elems, size_t type_size, size_t lane_count) const;

		/**
		 * Query whether profiling is enabled on this device.
		 */
		bool is_profiling_enabled() const noexcept;

		/**
		 * Make sure all operations on this device are finished.
		 */
		void synchronize() const;

		/**
		 * Get the profiling data for all executions of the given kernel on this device.
		 *
		 * \param kernel The kernel for which to get the profiling data
		 * \return Profiling data for the given kernel on this device
		 */
		ProfilingData getProfilingData(const cl_kernel& kernel) noexcept;

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
		 * Whether this device prefers SOA storage
		 *
		 * @todo This should be datatype specific
		 */
		const bool prefers_soa;

		/**
		 * The name of this device
		 */
		const std::string name;

		/**
		 * Whether profiling is enabled
		 */
		const bool profiling_enabled;

		/**
		 * Allow easy use of the command queue
		 */
		operator cl_command_queue() const noexcept;

		/**
		 * Allow easy use of the command queue
		 */
		cl_command_queue get_queue() const noexcept;

		/**
		 * Kernel profiling data
		 */
		std::map<cl_kernel, ProfilingData> profiling_data;
	};
}

#endif /* _HARDWARE_DEVICE_HPP_ */
