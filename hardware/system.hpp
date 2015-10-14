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
#include "hardwareParameters.hpp"
#include "../common_header_files/types.h"
#include "size_4.hpp"
#include <map>
#include <memory>
#include <tuple>

/**
 * This namespace contains potentially hardware specific code managing the system,
 * its devices and the physical representations of data as well as their actual manipulations.
 */
namespace hardware {

	/**
	 * @todo merge with other exceptions
	 */ 
	class OpenclException {
		public:
			OpenclException(int clerr);
			OpenclException(int clerr, std::string clname);
			OpenclException(int clerr, std::string clname, std::string filename, int linenumber);
			std::string what();
			int errorCode;
		private:
			std::string error_message;
	};

	class Device;
	class Transfer;
	class OpenClCode;

	class System {

	public:
		/**
		 * Create a new system representation.
		 * You should usually only do this once per application.
		 */
		explicit System(const meta::Inputparameters& parameters);
		System(const hardware::HardwareParametersInterface & systemParameters, const hardware::OpenClCode & kernelBuilder);

		~System();

		const std::vector<Device*>& get_devices() const noexcept;
		const meta::Inputparameters& get_inputparameters() const noexcept;
		const hardware::HardwareParametersInterface * getHardwareParameters() const noexcept;

		// non-copyable
		System& operator=(const System&) = delete;
		System(const System&) = delete;
		System() = delete;
		cl_context getContext() const;

		/**
		 * Get the size of the device grid.
		 */
		size_4 get_grid_size();

		Transfer * get_transfer(size_t from, size_t to, unsigned id) const;
		cl_platform_id get_platform() const;

	private:

		const meta::Inputparameters * params;
		std::vector<Device*> devices;
		cl_context context;
		cl_platform_id platform;
		size_4 grid_size;

		mutable std::map<std::tuple<size_t,size_t,unsigned>,std::unique_ptr<Transfer>> transfer_links;

		void initOpenCLPlatforms();
		void initOpenCLContext();
		void initOpenCLDevices();
		const hardware::HardwareParametersInterface * hardwareParameters;
		const hardware::OpenClCode * kernelBuilder;
		bool temporaryFlagForSystemConstructorVersion = false;
	};

	/**
	 * Print the profiling information of kernels on this system on any device.
	 */
	void print_profiling(const System& system, const std::string& filenameToWriteTo);
	void print_profiling(const System* system, const std::string& filenameToWriteTo);
}

#endif /* _HARDWARE_SYSTEM_HPP_ */
