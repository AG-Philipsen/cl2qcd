/** @file
 * Declaration of the hardware::System class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#ifndef _HARDWARE_SYSTEM_HPP_
#define _HARDWARE_SYSTEM_HPP_

#include "../meta/inputparameters.hpp"
#include "device.hpp"

/**
 * This namespace contains potentially hardware specific code managing the system,
 * it's devices and the physical representations of data as well as their actual manipulations.
 */
namespace hardware {

	/**
	 * OpenCL exception class,
	 * thrown if an opencl error occurs.
	 *
	 * @todo merge with other one
	 */ 
	class OpenclException {
		public:
			OpenclException(int clerr);
			OpenclException(int clerr, std::string clname);
			OpenclException(int clerr, std::string clname, std::string filename, int linenumber);
			std::string what();
		private:
			std::string error_message;
	};

	/**
	 * Representation of the system we are working on.
	 *
	 * Allows the querying and manipulation of the system and
	 * its hardware.
	 */
	class System {

	public:
		/**
		 * Create a new system representation.
		 * You should usually only do this once per application.
		 *
		 * \param inputparams The inputparameters of the application
		 * \param enableProfiling Whether we want to use profiling
		 */
		System(const meta::Inputparameters&, bool enable_profiling = false);

		~System();

		/**
		 * Get the devices of the system.
		 */
		const std::vector<Device*>& get_devices() const noexcept;

		/**
		 * Get the inputparameters of the system.
		 */
		const meta::Inputparameters& get_inputparameters() const noexcept;

		// non-copyable
		System& operator=(const System&) = delete;
		System(const System&) = delete;
		System() = delete;

		/**
		 * Allow to use the platform id
		 *
		 * \deprected This is not meant for wider application and only there to ease transition to the new architecture
		 */
		operator const cl_platform_id&() const noexcept;

		/**
		 * Allow to use the context
		 *
		 * \deprected This is not meant for wider application and only there to ease transition to the new architecture
		 */
		operator const cl_context&() const noexcept;

	private:

		/**
		 * The input paramters of the application
		 */
		const meta::Inputparameters& params;

		/**
		 * The devices of the system
		 */
		std::vector<Device*> devices;

		/**
		 * The OpenCL context used by the application
		 */
		cl_context context;

		/**
		 * The platform used by this system.
		 */
		cl_platform_id platform;
	};

}

#endif /* _HARDWARE_SYSTEM_HPP_ */
