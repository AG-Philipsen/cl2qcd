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
#include "openClKernelParameters.hpp"
#include "../common_header_files/types.h"
#include "../geometry/latticeGrid.hpp"
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
		System(const hardware::HardwareParametersInterface &, const hardware::code::OpenClKernelParametersInterface &);
		System(meta::Inputparameters&); //@todo: only for compatibility, remove!

		~System();

		const std::vector<Device*>& get_devices() const noexcept;
		const meta::Inputparameters& get_inputparameters() const noexcept; //@todo: remove
		const hardware::HardwareParametersInterface * getHardwareParameters() const noexcept;

		// non-copyable
		System& operator=(const System&) = delete;
		System(const System&) = delete;
		System() = delete;
		cl_context getContext() const;

		Transfer * get_transfer(size_t from, size_t to, unsigned id) const;
		cl_platform_id get_platform() const;

	private:
		std::vector<Device*> devices;
		cl_context context;
		cl_platform_id platform;
		LatticeGrid lG;

		mutable std::map<std::tuple<size_t,size_t,unsigned>,std::unique_ptr<Transfer>> transfer_links;

		void initOpenCLPlatforms();
		void initOpenCLContext();
		void initOpenCLDevices();
		const hardware::HardwareParametersInterface * hardwareParameters;
		const hardware::code::OpenClKernelParametersInterface * kernelParameters;
		const hardware::OpenClCode * kernelBuilder;
		//Remove when dependence on meta::InputParameters has been removed from physics package, for the moment it is needed in tests in physics!!!!
		const meta::Inputparameters& inputparameters;
	};

	/**
	 * Print the profiling information of kernels on this system on any device.
	 */
	void print_profiling(const System& system, const std::string& filenameToWriteTo);
	void print_profiling(const System* system, const std::string& filenameToWriteTo);
}

#endif /* _HARDWARE_SYSTEM_HPP_ */
