/** @file
 * Unit test for the physics::lattices::Gaugefield class
 */

#include "gaugefield.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::Gaugefield
#include <boost/test/unit_test.hpp>

#include "../../logger.hpp"
#include "../../meta/type_ops.hpp"
#include <stdexcept>

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	logger.debug() << "Devices: " << system.get_devices().size();
	physics::PRNG prng(system);

	// init from file
	Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
	BOOST_CHECK_CLOSE(gf.plaquette(), 1., 0.1);

	// init hot
	Gaugefield gf2(system, prng, true);

	// init cold
	Gaugefield gf3(system, prng, false);
	BOOST_CHECK_CLOSE(gf3.plaquette(), 1., 0.1);
}

BOOST_AUTO_TEST_CASE(save)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	Gaugefield gf(system, prng, true);
	gf.save("conf.test", 0);

	Gaugefield reread(system, prng, "conf.test");

	hmc_float orig_plaq, reread_plaq;
	hmc_float orig_tplaq, reread_tplaq;
	hmc_float orig_splaq, reread_splaq;
	hmc_complex orig_pol, reread_pol;

	gf.gaugeobservables(&orig_plaq, &orig_tplaq, &orig_splaq, &orig_pol);
	gf.gaugeobservables(&reread_plaq, &reread_tplaq, &reread_splaq, &reread_pol);

	BOOST_CHECK_EQUAL(orig_plaq, reread_plaq);
	BOOST_CHECK_EQUAL(orig_splaq, reread_splaq);
	BOOST_CHECK_EQUAL(orig_tplaq, reread_tplaq);
	BOOST_CHECK_EQUAL(orig_pol, reread_pol);
}

BOOST_AUTO_TEST_CASE(rectangles)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo", "--gaugeact=wilson"};
	meta::Inputparameters params(2, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
	BOOST_CHECK_THROW(gf.rectangles(), std::logic_error);

	const char * _params2[] = {"foo", "--gaugeact=tlsym"};
	meta::Inputparameters params2(2, _params2);
	hardware::System system2(params2);
	physics::PRNG prng2(system2);

	Gaugefield gf2(system2, prng2, std::string(SOURCEDIR) + "/tests/conf.00200");
	BOOST_CHECK_CLOSE(gf2.rectangles(), 6144, 0.1);
}
