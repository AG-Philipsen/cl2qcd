/** @file
 * Unit test for the physics::lattices::Spinorfield_eo class
 */

#include "spinorfield_eo.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::Spinorfield_eo
#include <boost/test/unit_test.hpp>

#include "../../logger.hpp"

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	logger.debug() << "Devices: " << system.get_devices().size();

	Spinorfield_eo sf(system);

	// init from file
	BOOST_ERROR("not implemented");
}

BOOST_AUTO_TEST_CASE(save)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	Spinorfield_eo sf(system);
	BOOST_ERROR("not implemented");
}
