/** @file
 * Declaration of the hardware::System class
 *
 * Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
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

#ifndef _HARDWARE_SYSTEM_HPP_
#define _HARDWARE_SYSTEM_HPP_

#include "../meta/inputparameters.hpp"
#include "../types.h"
#include "../meta/size_4.hpp"
#include <map>
#include <memory>
#include <tuple>

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

	// forward decleartion of used class
	class Device;
	class Transfer;

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
		explicit System(const meta::Inputparameters&, bool enable_profiling = false);

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
		 * Allow to use the context
		 *
		 * \deprected This is not meant for wider application and only there to ease transition to the new architecture
		 */
		operator const cl_context&() const noexcept;

		/**
		 * Get the size of the device grid.
		 */
		size_4 get_grid_size();

		Transfer * get_transfer(size_t from, size_t to, unsigned id) const;

		/**
		 * Get the platform used.
		 */
		cl_platform_id get_platform() const;

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

		/**
		 * The size of the device grid.
		 */
		size_4 grid_size;

		mutable std::map<std::tuple<size_t,size_t,unsigned>,std::unique_ptr<Transfer>> transfer_links;
	};

	/**
	 * Print the profiling information of kernels on this system on any device.
	 *
	 * \param system The system of which to print the profiling info
	 * \param filename The file to write the profiling information to
	 */
	void print_profiling(const System& system, const std::string& filename);

	/**
	 * Print the profiling information of kernels on this system on any device.
	 *
	 * \param system The system of which to print the profiling info
	 * \param filename The file to write the profiling information to
	 */
	void print_profiling(const System* system, const std::string& filename);
}

#endif /* _HARDWARE_SYSTEM_HPP_ */
