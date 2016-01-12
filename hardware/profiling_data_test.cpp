/** @file
 * Testcases for the hardware::ProfilingData class
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
#define BOOST_TEST_MODULE hardware::ProfilingData
#include <boost/test/unit_test.hpp>

#include "profiling_data.hpp"

#include "system.hpp"
#include "device.hpp"
#include "code/buffer.hpp"
#include "interfaceMockups.hpp"

BOOST_AUTO_TEST_CASE(initial)
{
	using namespace hardware;

	const ProfilingData data;
	BOOST_CHECK_EQUAL(data.get_total_time(), 0);
	BOOST_CHECK_EQUAL(data.get_num_values(), 0);
}

namespace hardware {
	cl_command_queue profilingDataTestCommandQueueHelper(const hardware::Device * device)
	{
		return device->get_queue();
	}
}

BOOST_AUTO_TEST_CASE(add_value)
{
	using namespace hardware;

	const hardware::HardwareParametersMockupWithProfiling hardwareParameters(4,4);
	const hardware::code::OpenClKernelParametersMockup kernelParameters(4,4);
	hardware::System system( hardwareParameters, kernelParameters );

	// there should always be at least one device
	// otherwise code or system is broken
	// in both cases it is good to fail
	const std::vector<Device*>& devices = system.get_devices();
	BOOST_REQUIRE_GE(devices.size(), 1);

	// query some data
	for(const Device * device : devices)
	{
		ProfilingData data;
		ProfilingData old_data;
		cl_event event;
		for(size_t i = 0; i < 3; i++)
		{
			clEnqueueMarker(profilingDataTestCommandQueueHelper(device), &event);
			clWaitForEvents(1, &event);
			data += event;
			clReleaseEvent(event);
			BOOST_REQUIRE_EQUAL(data.get_num_values(), old_data.get_num_values() + 1);
			BOOST_REQUIRE_GE(data.get_total_time(), old_data.get_total_time());
			old_data = data;
		}
	}
}
