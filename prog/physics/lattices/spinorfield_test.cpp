/** @file
 * Unit test for the physics::lattices::Spinorfield class
 */

#include "spinorfield.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::Spinorfield
#include <boost/test/unit_test.hpp>

#include "../../logger.hpp"

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	logger.debug() << "Devices: " << system.get_devices().size();

	Spinorfield sf(system);

	// init from file
	BOOST_ERROR("not implemented");
}

BOOST_AUTO_TEST_CASE(save)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);

	Spinorfield sf(system);
	BOOST_ERROR("not implemented");
}

BOOST_AUTO_TEST_CASE(fields)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);

	auto fields = create_spinorfields(system, 2);

	release_spinorfields(fields);
}
