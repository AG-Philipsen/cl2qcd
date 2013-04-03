/** @file
 * Tests of the fermion force algorithms
 */

#include "fermion_force.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::fermion_force
#include <boost/test/unit_test.hpp>

#include "../lattices/util.hpp"

BOOST_AUTO_TEST_CASE(fermion_force)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=16"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);

		Gaugefield gf(system, prng, false);
		Spinorfield sf1(system);
		Spinorfield sf2(system);
		Gaugemomenta gm(system);

		pseudo_randomize<Spinorfield, spinor>(&sf1, 11);
		pseudo_randomize<Spinorfield, spinor>(&sf2, 12);
		gm.zero();

		physics::algorithms::fermion_force(&gm, sf1, sf2, gf);
		BOOST_CHECK_CLOSE(squarenorm(gm), 18550.897680606064, 0.01);
	}

	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);

		Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
		Spinorfield sf1(system);
		Spinorfield sf2(system);
		Gaugemomenta gm(system);

		pseudo_randomize<Spinorfield, spinor>(&sf1, 13);
		pseudo_randomize<Spinorfield, spinor>(&sf2, 14);
		gm.zero();

		physics::algorithms::fermion_force(&gm, sf1, sf2, gf);
		BOOST_CHECK_CLOSE(squarenorm(gm), 3561.5546616746424, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(fermion_force_eo)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=16"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);

		Gaugefield gf(system, prng, false);
		Spinorfield src(system);
		Spinorfield_eo sf1(system);
		Spinorfield_eo sf2(system);
		Gaugemomenta gm(system);

		pseudo_randomize<Spinorfield, spinor>(&src, 15);
		convert_to_eoprec(&sf1, &sf2, src);
		gm.zero();

		physics::algorithms::fermion_force(&gm, sf1, sf2, EVEN, gf);
		BOOST_CHECK_CLOSE(squarenorm(gm), 6180.0548464319563, 0.01);
		physics::algorithms::fermion_force(&gm, sf1, sf2, ODD, gf);
		BOOST_CHECK_CLOSE(squarenorm(gm), 18590.240819832132, 0.01);
	}

	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);

		Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
		Spinorfield src(system);
		Spinorfield_eo sf1(system);
		Spinorfield_eo sf2(system);
		Gaugemomenta gm(system);

		pseudo_randomize<Spinorfield, spinor>(&src, 16);
		convert_to_eoprec(&sf1, &sf2, src);
		gm.zero();

		physics::algorithms::fermion_force(&gm, sf1, sf2, EVEN, gf);
		BOOST_CHECK_CLOSE(squarenorm(gm), 1294.880037707632, 0.01);
		physics::algorithms::fermion_force(&gm, sf1, sf2, ODD, gf);
		BOOST_CHECK_CLOSE(squarenorm(gm), 3659.59932413153, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(calc_fermion_forces)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);

		Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
		Spinorfield sf1(system);
		Gaugemomenta gm(system);

		pseudo_randomize<Spinorfield, spinor>(&sf1, 21);
		gm.zero();

		physics::algorithms::calc_fermion_forces(&gm, gf, sf1, system);
		BOOST_CHECK_CLOSE(squarenorm(gm), 42199.514415107173, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(calc_fermion_forces_eo)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);

		Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
		Spinorfield src(system);
		Spinorfield_eo sf1(system);
		Spinorfield_eo sf2(system);
		Gaugemomenta gm(system);

		pseudo_randomize<Spinorfield, spinor>(&src, 22);
		convert_to_eoprec(&sf1, &sf2, src);
		gm.zero();

		physics::algorithms::calc_fermion_forces(&gm, gf, sf1, system);
		BOOST_CHECK_CLOSE(squarenorm(gm), 3441.344988280136, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(calc_detratio_forces)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);

		Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
		Spinorfield sf1(system);
		Gaugemomenta gm(system);

		pseudo_randomize<Spinorfield, spinor>(&sf1, 21);
		gm.zero();

		physics::algorithms::calc_detratio_forces(&gm, gf, sf1, system);
		BOOST_CHECK_CLOSE(squarenorm(gm), 2.8236650583738432e-12, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(calc_detratio_forces_eo)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);

		Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
		Spinorfield src(system);
		Spinorfield_eo sf1(system);
		Spinorfield_eo sf2(system);
		Gaugemomenta gm(system);

		pseudo_randomize<Spinorfield, spinor>(&src, 22);
		convert_to_eoprec(&sf1, &sf2, src);
		gm.zero();

		physics::algorithms::calc_detratio_forces(&gm, gf, sf1, system);
		BOOST_CHECK_CLOSE(squarenorm(gm), 3.2577363488458202e-12, 0.01);
	}
}
