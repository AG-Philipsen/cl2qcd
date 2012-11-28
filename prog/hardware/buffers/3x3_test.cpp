/** @file
 * Testcases for the hardware::buffers::Matrix3x3 class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "3x3.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hardware::buffers::Matrix3x3
#include <boost/test/unit_test.hpp>

#include "../system.hpp"
#include "../../meta/util.hpp"
#include "../../meta/type_ops.hpp"

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace hardware;

	System system(meta::Inputparameters(0, 0));
for(Device * device : system.get_devices()) {

		hardware::buffers::Matrix3x3 dummy(meta::get_vol4d(system.get_inputparameters()) * NDIM, device);
		const cl_mem * tmp = dummy;
		BOOST_CHECK(tmp);
		BOOST_CHECK(*tmp);
	}
}
