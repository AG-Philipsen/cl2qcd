/** @file
 * Testcases for the hardware::buffers::Spinor class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "su3vec.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hardware::buffers::SU3vec
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

		SU3vec dummy(meta::get_vol4d(system.get_inputparameters()), device);
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
	const size_t elems = meta::get_vol4d(system.get_inputparameters()) / 2;
for(Device * device : system.get_devices()) {
		su3vec* buf = new su3vec[elems];
		su3vec* buf2 = new su3vec[elems];
		SU3vec dummy(elems, device);
		fill(buf, elems, 1);
		fill(buf2, elems, 2);
		dummy.load(buf);
		dummy.dump(buf2);
		BOOST_CHECK_EQUAL_COLLECTIONS(buf, buf + elems, buf2, buf2 + elems);
		delete[] buf;
		delete[] buf2;
	}
}

BOOST_AUTO_TEST_CASE(copy)
{
	using namespace hardware;
	using namespace hardware::buffers;

	System system(meta::Inputparameters(0, 0));
	const size_t elems = meta::get_vol4d(system.get_inputparameters()) / 2;
for(Device * device : system.get_devices()) {
		su3vec* buf = new su3vec[elems];
		su3vec* buf2 = new su3vec[elems];
		SU3vec dummy(elems, device);
		SU3vec dummy2(elems, device);

		fill(buf, elems, 1);
		fill(buf2, elems, 2);
		dummy.load(buf);
		copyData(&dummy2, &dummy);
		dummy2.dump(buf2);
		BOOST_CHECK_EQUAL_COLLECTIONS(buf, buf + elems, buf2, buf2 + elems);

		delete[] buf;
		delete[] buf2;
	}
}
