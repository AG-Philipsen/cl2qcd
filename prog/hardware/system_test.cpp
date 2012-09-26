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

	BOOST_REQUIRE(static_cast<const cl_platform_id>(system));
	BOOST_REQUIRE(static_cast<const cl_context>(system));

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
	meta::Inputparameters params3(2, _params3);
	System system3(params3);
	const std::vector<Device*>& devices3 = system3.get_devices();
	BOOST_REQUIRE_EQUAL(devices3.size(), 1);

	// check that only device 0 is returned
	const char * _params4[] = {"foo", "--device=0"};
	meta::Inputparameters params4(2, _params4);
	System system4(params4);
	const std::vector<Device*>& devices4 = system4.get_devices();
	BOOST_REQUIRE_EQUAL(devices4.size(), 1);

	// check whether GPUs/CPUs can be disabled
	const char * _params5[] = {"foo", "--use_gpu=false"};
	meta::Inputparameters params5(2, _params5);
	System system5(params5);
	const std::vector<Device*>& devices5 = system5.get_devices();
for(Device * device : devices5) {
		BOOST_REQUIRE_NE(device->get_device_type(), CL_DEVICE_TYPE_GPU);
	}
	const char * _params6[] = {"foo", "--use_cpu=false"};
	meta::Inputparameters params6(2, _params6);
	System system6(params6);
	const std::vector<Device*>& devices6 = system6.get_devices();
for(Device * device : devices6) {
		BOOST_REQUIRE_NE(device->get_device_type(), CL_DEVICE_TYPE_CPU);
	}
}

BOOST_AUTO_TEST_CASE(dump_source_if_debugging)
{
	BOOST_FAIL("not implemented");
}
