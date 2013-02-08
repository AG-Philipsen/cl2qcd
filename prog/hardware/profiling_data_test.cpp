/** @file
 * Testcases for the hardware::ProfilingData class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hardware::ProfilingData
#include <boost/test/unit_test.hpp>

#include "profiling_data.hpp"

#include "system.hpp"
#include "device.hpp"
#include "code/buffer.hpp"

BOOST_AUTO_TEST_CASE(initial)
{
	using namespace hardware;

	const ProfilingData data;
	BOOST_CHECK_EQUAL(data.get_total_time(), 0);
	BOOST_CHECK_EQUAL(data.get_num_values(), 0);
}

namespace hardware {
cl_command_queue profiling_data_test_command_queue_helper(const hardware::Device * device)
{
	return device->get_queue();
}
}

BOOST_AUTO_TEST_CASE(add_value)
{
	using namespace hardware;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	System system(params, true);

	// there should always be at least one device
	// otherwise code or system is broken
	// in both cases it is good to fail
	const std::vector<Device*>& devices = system.get_devices();
	BOOST_REQUIRE_GE(devices.size(), 1);

	// query some data
for(const Device * device : devices) {

		ProfilingData data;
		ProfilingData old_data;
		cl_event event;
		for(size_t i = 0; i < 3; i++) {
			clEnqueueMarker(profiling_data_test_command_queue_helper(device), &event);
			clWaitForEvents(1, &event);
			data += event;
			clReleaseEvent(event);
			BOOST_REQUIRE_EQUAL(data.get_num_values(), old_data.get_num_values() + 1);
			BOOST_REQUIRE_GE(data.get_total_time(), old_data.get_total_time());
			old_data = data;
		}
	}
}
