/** @file
 * Testcases for the hardware::buffers::SU3 class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "su3.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hardware::buffers::SU3
#include <boost/test/unit_test.hpp>

#include "../system.hpp"
#include "../../meta/util.hpp"

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace hardware;
	using namespace hardware::buffers;

	System system(meta::Inputparameters(0, 0));
for(Device * device : system.get_devices()) {

		SU3 dummy(meta::get_vol4d(system.get_inputparameters()), device);
		const cl_mem * tmp = dummy;
		BOOST_CHECK(tmp);
		BOOST_CHECK(*tmp);
	}
}

BOOST_AUTO_TEST_CASE(import_export)
{
	BOOST_FAIL("not implemented");
}
