/** @file
 * Testcases for the hardware::buffers::Buffer class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "prng_buffer.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hardware::buffers::PRNGBuffer
#include <boost/test/unit_test.hpp>

#include "../system.hpp"

BOOST_AUTO_TEST_CASE(get_prng_buffer_size)
{
	using namespace hardware;
	using namespace hardware::buffers;

	System system(meta::Inputparameters(0, 0));
for(Device * device : system.get_devices()) {

		int elems = hardware::buffers::get_prng_buffer_size(device);

		BOOST_CHECK_GT(elems, 0);
		BOOST_CHECK_LT(elems, 1e6);
	}
}

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace hardware;
	using namespace hardware::buffers;

	System system(meta::Inputparameters(0, 0));
for(Device * device : system.get_devices()) {

		PRNGBuffer dummy(device);
		BOOST_CHECK_EQUAL(dummy.get_elements(), hardware::buffers::get_prng_buffer_size(device));
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
for(Device * device : system.get_devices()) {

		PRNGBuffer buffer(device);
		int elems = buffer.get_elements();
		PRNGBuffer::prng_state_t * in = new PRNGBuffer::prng_state_t[elems];
		// TODO fill with random data
		buffer.load(in);
		PRNGBuffer::prng_state_t * out = new PRNGBuffer::prng_state_t[elems];
		buffer.dump(out);

		BOOST_CHECK_EQUAL_COLLECTIONS(reinterpret_cast<uint64_t*>(in), reinterpret_cast<uint64_t*>(in + elems), reinterpret_cast<uint64_t*>(out), reinterpret_cast<uint64_t*>(out + elems));

		delete[] out;
		delete[] in;
	}
}
