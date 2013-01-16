/** @file
 * Unit test for the physics::lattices::Scalar<hmc_complex> class
 */

#include "scalar_complex.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::Scalar<hmc_complex>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>

#include "../../logger.hpp"
#include "../../meta/type_ops.hpp"

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	logger.debug() << "Devices: " << system.get_devices().size();

	Scalar<hmc_complex> foo(system);
	hmc_complex ref = {42., 13.};
	foo.store(ref);

	BOOST_REQUIRE_EQUAL(ref, foo.get());
}

BOOST_AUTO_TEST_CASE(multiplication)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	logger.debug() << "Devices: " << system.get_devices().size();

	Scalar<hmc_complex> left(system);
	Scalar<hmc_complex> right(system);
	Scalar<hmc_complex> res(system);

	left.store( {1., 2.});
	right.store( {3., 4.});

	multiply(&res, left, right);

	hmc_complex ref = { -5., 10.};
	BOOST_REQUIRE_EQUAL(ref, res.get());
}

BOOST_AUTO_TEST_CASE(division)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	logger.debug() << "Devices: " << system.get_devices().size();

	Scalar<hmc_complex> left(system);
	Scalar<hmc_complex> right(system);
	Scalar<hmc_complex> res(system);

	left.store( {1., 2.});
	right.store( {3., 4.});

	divide(&res, left, right);

	hmc_complex res_host = res.get();
	BOOST_REQUIRE_CLOSE(0.44, res_host.re, 1);
	BOOST_REQUIRE_CLOSE(0.08, res_host.im, 1);
}
