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

	{
		const char * _params[] = {"foo", "--ntime=8"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		logger.debug() << "Devices: " << system.get_devices().size();
		physics::PRNG prng(system);

		// init hot
		Gaugefield gf2(system, prng, true);

		// init cold
		Gaugefield gf3(system, prng, false);
		BOOST_CHECK_CLOSE(gf3.plaquette(), 1., 0.1);
	}

	{
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		logger.debug() << "Devices: " << system.get_devices().size();
		physics::PRNG prng(system);

		// init from file
		Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
		BOOST_CHECK_CLOSE(gf.plaquette(), 0.57107711169452713, 0.1);
	}
}

BOOST_AUTO_TEST_CASE(save)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo"};
	meta::Inputparameters params(1, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	Gaugefield gf(system, prng, false);
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

	const char * _params[] = {"foo", "--gaugeact=wilson", "--ntime=4"};
	meta::Inputparameters params(3, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
	BOOST_CHECK_THROW(gf.rectangles(), std::logic_error);

	const char * _params2[] = {"foo", "--gaugeact=tlsym", "--ntime=4"};
	meta::Inputparameters params2(3, _params2);
	hardware::System system2(params2);
	physics::PRNG prng2(system2);

	Gaugefield gf2(system2, prng2, std::string(SOURCEDIR) + "/tests/conf.00200");
	BOOST_CHECK_CLOSE(gf2.rectangles(), 1103.2398401620451, 0.1);
}

BOOST_AUTO_TEST_CASE(polyakov)
{
	using namespace physics::lattices;

	const char * _params[] = {"foo", "--ntime=4"};
	meta::Inputparameters params(2, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
	hmc_complex pol = gf.polyakov();
	BOOST_CHECK_CLOSE(pol.re, -0.11349672123636857, 0.1);
	BOOST_CHECK_CLOSE(pol.im, 0.22828243566855227, 0.1);
}

BOOST_AUTO_TEST_CASE(halo_update)
{
	using namespace physics::lattices;

	hmc_float orig_plaq, new_plaq;
	hmc_float orig_tplaq, new_tplaq;
	hmc_float orig_splaq, new_splaq;
	hmc_complex orig_pol, new_pol;

	// simple test, gaugeobservables should not get changed by halo exchange
	// if the original gaugefield is given
	{
		const char * _params[] = {"foo", "--ntime=16"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);

		Gaugefield gf(system, prng, false);
		gf.gaugeobservables(&orig_plaq, &orig_tplaq, &orig_splaq, &orig_pol);
		gf.update_halo();
		gf.gaugeobservables(&new_plaq, &new_tplaq, &new_splaq, &new_pol);

		BOOST_CHECK_EQUAL(orig_plaq, new_plaq);
		BOOST_CHECK_EQUAL(orig_splaq, new_splaq);
		BOOST_CHECK_EQUAL(orig_tplaq, new_tplaq);
		BOOST_CHECK_EQUAL(orig_pol, new_pol);
	}

	// in a hot initialization halo cells will have wrong values until they are exchanged,
	// but only if there is more than one device
	// TODO hot initalization should automatically make sure the halo get's updated
	{
		const char * _params[] = {"foo"};
		meta::Inputparameters params(1, _params);
		hardware::System system(params);
		physics::PRNG prng(system);

		Gaugefield gf(system, prng, true);
		gf.gaugeobservables(&orig_plaq, &orig_tplaq, &orig_splaq, &orig_pol);
		gf.update_halo();
		gf.gaugeobservables(&new_plaq, &new_tplaq, &new_splaq, &new_pol);

		if(system.get_devices().size() > 1) {
			BOOST_CHECK_NE(orig_plaq, new_plaq);
			BOOST_CHECK_EQUAL(orig_splaq, new_splaq);
			BOOST_CHECK_NE(orig_tplaq, new_tplaq);
		} else {
			BOOST_CHECK_EQUAL(orig_plaq, new_plaq);
			BOOST_CHECK_EQUAL(orig_splaq, new_splaq);
			BOOST_CHECK_EQUAL(orig_tplaq, new_tplaq);
			BOOST_CHECK_EQUAL(orig_pol, new_pol);
		}
	}

	{
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);

		Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
		gf.gaugeobservables(&orig_plaq, &orig_tplaq, &orig_splaq, &orig_pol);
		gf.update_halo();
		gf.gaugeobservables(&new_plaq, &new_tplaq, &new_splaq, &new_pol);

		BOOST_CHECK_EQUAL(orig_plaq, new_plaq);
		BOOST_CHECK_EQUAL(orig_splaq, new_splaq);
		BOOST_CHECK_EQUAL(orig_tplaq, new_tplaq);
		BOOST_CHECK_EQUAL(orig_pol, new_pol);
	}

}
