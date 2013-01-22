/** @file
 * Unit test for the physics::lattices::Gaugemomenta class
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "gaugemomenta.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::Gaugemomenta
#include <boost/test/unit_test.hpp>

#include "../../logger.hpp"

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	logger.debug() << "Devices: " << system.get_devices().size();

	Gaugemomenta gf(system);

	BOOST_REQUIRE_NE(gf.get_buffers().size(), 0u);
}
