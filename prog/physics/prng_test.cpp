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
