/** @file
 * Unit test for the physics::lattices::Gaugefield class
 */

#include "gaugefield.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::Gaugefield
#include <boost/test/unit_test.hpp>

#include "../../logger.hpp"

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	logger.debug() << "Devices: " << system.get_devices().size();
	physics::PRNG prng(system);

	// init from file
	BOOST_ERROR("not implemented");

	// init hot
	Gaugefield gf2(system, prng, true);
	// TODO test
	// init cold
	Gaugefield gf3(system, prng);
	// TODO test
}

BOOST_AUTO_TEST_CASE(save)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	Gaugefield gf(system, prng);
	BOOST_ERROR("not implemented");
}
