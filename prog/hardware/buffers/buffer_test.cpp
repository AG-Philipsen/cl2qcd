/** @file
 * Testcases for the hardware::buffers::Buffer class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "buffer.hpp"
#include "../system.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hardware::buffers::Buffer
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace hardware;
	using namespace hardware::buffers;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	System system(params);
	const std::vector<Device*>& devices = system.get_devices();
for(Device * device : devices) {

		Buffer dummy(sizeof(float), device);
		BOOST_REQUIRE_EQUAL(dummy.get_bytes(), sizeof(float));
		const cl_mem * tmp = dummy;
		BOOST_REQUIRE(tmp);
		BOOST_REQUIRE(*tmp);
	}
}
