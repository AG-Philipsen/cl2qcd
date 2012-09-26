/** @file
 * Testcases for the hardware::Device class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hardware::Device
#include <boost/test/unit_test.hpp>

#include "system.hpp"

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
for(Device * device : devices) {
		device->is_double_supported();
	}
}

BOOST_AUTO_TEST_CASE(compile)
{
	BOOST_FAIL("not implemented");
}

BOOST_AUTO_TEST_CASE(profiling)
{
	BOOST_FAIL("not implemented");
}
