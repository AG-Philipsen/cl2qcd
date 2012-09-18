/** @file
 * Testcases for the hardware::System class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hardware::System
#include <boost/test/unit_test.hpp>

#include "system.hpp"

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace hardware;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);

	// test it doesn't blow up
	System system(params);

	// and the inputparameters should be the provided ones
	BOOST_REQUIRE_EQUAL(&system.get_inputparameters(), &params);
}

BOOST_AUTO_TEST_CASE(devices)
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

	// if double is required, check that all devices support double
for(Device * device : devices) {
		BOOST_REQUIRE_EQUAL(device->is_double_supported(), true);
	}
//	const char * _params2[] = {"foo"};
//	meta::Inputparameters params2(2, _params2);
//	BOOST_FAIL("not implemented");

	// check that nor more than the specified amount of devices is returned
	const char * _params3[] = {"foo", "--num_dev=1"};
	meta::Inputparameters params3(1, _params3);
	System system3(params3);
	const std::vector<Device*>& devices3 = system3.get_devices();
	BOOST_REQUIRE_EQUAL(devices3.size(), 1);

	// check that only device 0 is returned
	const char * _params4[] = {"foo", "--device=0"};
	meta::Inputparameters params4(1, _params4);
	System system4(params4);
	const std::vector<Device*>& devices4 = system4.get_devices();
	BOOST_REQUIRE_EQUAL(devices4.size(), 1);
}

BOOST_AUTO_TEST_CASE(dump_source_if_debugging)
{
	BOOST_FAIL("not implemented");
}
