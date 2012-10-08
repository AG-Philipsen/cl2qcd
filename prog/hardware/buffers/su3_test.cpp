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
#include "../../meta/type_ops.hpp"

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace hardware;
	using namespace hardware::buffers;

	System system(meta::Inputparameters(0, 0));
for(Device * device : system.get_devices()) {

		SU3 dummy(meta::get_vol4d(system.get_inputparameters()) * NDIM, device);
		const cl_mem * tmp = dummy;
		BOOST_CHECK(tmp);
		BOOST_CHECK(*tmp);
	}
}

BOOST_AUTO_TEST_CASE(import_export)
{
	using namespace hardware;
	using namespace hardware::buffers;

	System system(meta::Inputparameters(0, 0));
	const size_t elems = meta::get_vol4d(system.get_inputparameters()) * NDIM;
for(Device * device : system.get_devices()) {
		Matrixsu3* buf(new Matrixsu3[elems]);
		Matrixsu3* buf2(new Matrixsu3[elems]);
		SU3 dummy(elems, device);
		if(dummy.is_soa()) {
			BOOST_CHECK_THROW(dummy.load(buf), std::logic_error);
			BOOST_CHECK_THROW(dummy.dump(buf), std::logic_error);
		} else {
			fill(buf, elems, 1);
			fill(buf, elems, 2);
			dummy.load(buf);
			dummy.dump(buf2);
			BOOST_CHECK_EQUAL_COLLECTIONS(buf, buf + elems, buf2, buf2 + elems);
		}
		delete[] buf;
		delete[] buf2;
	}
}
