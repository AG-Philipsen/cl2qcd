/** @file
 * Tests of the physics::lattices::SwappableSpinorfield class
 */

#include "swappable_spinorfield.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::SwappableSpinorfield
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	logger.debug() << "Devices: " << system.get_devices().size();

	SwappableSpinorfield sf(system);
}

BOOST_AUTO_TEST_CASE(swap)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	logger.debug() << "Devices: " << system.get_devices().size();

	SwappableSpinorfield sf(system);
	BOOST_REQUIRE_LT(0, sf.get_buffers().size());

	sf.swap_out();
	BOOST_REQUIRE_EQUAL(0, sf.get_buffers().size());

	sf.swap_in();
	BOOST_REQUIRE_LT(0, sf.get_buffers().size());
}

