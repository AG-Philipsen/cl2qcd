/** @file
 * Testcases for the hardware::Device class
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

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hardware::Device
#include <boost/test/unit_test.hpp>

#include "system.hpp"
#include "device.hpp"

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace hardware;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	System system(params);

	// there should always be at least one device
	// otherwise code or system is broken
	// in both cases it is good to fail
	const std::vector<Device*>& devices = system.get_devices();
	BOOST_REQUIRE_GE(devices.size(), 1);

	// query some data
for(const Device * device : devices) {
		device->is_double_supported();
		BOOST_REQUIRE_NE(device->get_preferred_local_thread_num(), 0);
		BOOST_REQUIRE_NE(device->get_preferred_global_thread_num(), 0);
		const size_t recommended_stride = device->recommend_stride(1024, 16, 2);
		BOOST_REQUIRE_GE(recommended_stride, 1024);
	}
}

BOOST_AUTO_TEST_CASE(compile)
{
	using namespace hardware;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	System system(params);

	// there should always be at least one device
	// otherwise code or system is broken
	// in both cases it is good to fail
	const std::vector<Device*>& devices = system.get_devices();
	BOOST_REQUIRE_GE(devices.size(), 1);

	// query some data
for(const Device * device : devices) {
		cl_kernel foo = device->create_kernel("foo", "") << "hardware/device_test.cl";
	}
}
