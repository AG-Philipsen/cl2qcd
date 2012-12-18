/** @file
 * Ranlux PRNG unit test
 */

#include "prng.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::PRNG
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace physics;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);

	PRNG prng(system);
	// TODO test initial state
}

BOOST_AUTO_TEST_CASE(use_on_devices)
{
	BOOST_FAIL("not implemented");
}

BOOST_AUTO_TEST_CASE(store_and_resume)
{
	using namespace physics;

	const char * _params[] = {"foo", "--host_seed=5"};
	meta::Inputparameters params(2, _params);
	hardware::System system(params);

	PRNG prng(system);
	prng.store("tmp.prngstate");

	double prng1_res = prng.get_double();

	const char * _params2[] = {"foo", "--initial_prng_state=tmp.prngstate"};
	meta::Inputparameters params2(2, _params2);
	hardware::System system2(params2);

	PRNG prng2(system2);

	double prng2_res = prng2.get_double();

	BOOST_CHECK_EQUAL(prng1_res, prng2_res);

	logger.info() << "Now checking buffers";

	// TODO test also for device
	for(size_t i = 0; i < prng.get_buffers().size(); ++i) {
		auto buf1 = prng.get_buffers().at(i);
		auto buf2 = prng2.get_buffers().at(i);

		char* prng1_state = new char[buf1->get_bytes()];
		buf1->dump(reinterpret_cast<hardware::buffers::PRNGBuffer::prng_state_t *>(prng1_state));

		char* prng2_state = new char[buf2->get_bytes()];
		buf2->dump(reinterpret_cast<hardware::buffers::PRNGBuffer::prng_state_t *>(prng2_state));

		BOOST_CHECK_EQUAL_COLLECTIONS(prng1_state, prng1_state + buf1->get_bytes(), prng2_state, prng2_state + buf2->get_bytes());
		logger.info() << "Checked buffer " << i;
	}
}
