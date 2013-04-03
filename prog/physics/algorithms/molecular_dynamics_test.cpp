/** @file
 * Tests of the molecular dynamics algorithms
 */

#include "molecular_dynamics.hpp"
#include "forces.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::molecular_dynamics
#include <boost/test/unit_test.hpp>

#include "../lattices/util.hpp"

BOOST_AUTO_TEST_CASE(md_update_gaugefield)
{
	using namespace physics::lattices;
	using namespace physics::algorithms;
	{
		const char * _params[] = {"foo", "--ntime=16"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);

		Gaugefield gf(system, prng, false);
		Gaugemomenta gm(system);
		gm.zero();

		hmc_float ref = gf.plaquette();
		gauge_force(&gm, gf);
		physics::algorithms::md_update_gaugefield(&gf, gm, .5);
		BOOST_CHECK_CLOSE(gf.plaquette(), ref, 0.01);
	}

	{
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);

		Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
		Gaugemomenta gm(system);
		gm.zero();

		gauge_force(&gm, gf);
		physics::algorithms::md_update_gaugefield(&gf, gm, .5);
		BOOST_CHECK_CLOSE(gf.plaquette(), 0.0060440132434446334, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(md_update_spinorfield)
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

		pseudo_randomize<Spinorfield, spinor>(&sf1, 11);
		pseudo_randomize<Spinorfield, spinor>(&sf2, 12);

		physics::algorithms::md_update_spinorfield(&sf2, gf, sf1, system);
		BOOST_CHECK_CLOSE(squarenorm(sf2), 2523.6197747176057, 0.01);
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

		physics::algorithms::md_update_spinorfield(&sf2, gf, sf1, system);
		BOOST_CHECK_CLOSE(squarenorm(sf2), 2588.2852881545778, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(md_update_spinorfield_eo)
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

		pseudo_randomize<Spinorfield, spinor>(&src, 15);
		convert_to_eoprec(&sf1, &sf2, src);

		physics::algorithms::md_update_spinorfield(&sf2, gf, sf1, system);
		BOOST_CHECK_CLOSE(squarenorm(sf2), 1122.194230730885, 0.01);
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

		pseudo_randomize<Spinorfield, spinor>(&src, 16);
		convert_to_eoprec(&sf1, &sf2, src);

		physics::algorithms::md_update_spinorfield(&sf2, gf, sf1, system);
		BOOST_CHECK_CLOSE(squarenorm(sf2), 1114.3019247079062, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(md_update_spinorfield_mp)
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

		pseudo_randomize<Spinorfield, spinor>(&sf1, 21);
		pseudo_randomize<Spinorfield, spinor>(&sf2, 22);

		physics::algorithms::md_update_spinorfield(&sf2, gf, sf1, system);
		BOOST_CHECK_CLOSE(squarenorm(sf2), 2539.2078579177487, 0.01);
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

		pseudo_randomize<Spinorfield, spinor>(&sf1, 23);
		pseudo_randomize<Spinorfield, spinor>(&sf2, 24);
		gm.zero();

		physics::algorithms::md_update_spinorfield(&sf2, gf, sf1, system);
		BOOST_CHECK_CLOSE(squarenorm(sf2), 2658.5475465525888, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(md_update_spinorfield_mp_eo)
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

		pseudo_randomize<Spinorfield, spinor>(&src, 25);
		convert_to_eoprec(&sf1, &sf2, src);

		physics::algorithms::md_update_spinorfield(&sf2, gf, sf1, system);
		BOOST_CHECK_CLOSE(squarenorm(sf2), 1105.9891302664669, 0.01);
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

		pseudo_randomize<Spinorfield, spinor>(&src, 26);
		convert_to_eoprec(&sf1, &sf2, src);

		physics::algorithms::md_update_spinorfield(&sf2, gf, sf1, system);
		BOOST_CHECK_CLOSE(squarenorm(sf2), 1081.4370785519473, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(md_update_gaugemomentum_gauge)
{
	using namespace physics::lattices;
	using namespace physics::algorithms;
	{
		const char * _params[] = {"foo", "--ntime=16"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);

		Gaugefield gf(system, prng, false);
		Gaugemomenta gm(system);
		gm.zero();

		physics::algorithms::md_update_gaugemomentum_gauge(&gm, .5, gf, system);
		BOOST_CHECK_CLOSE(squarenorm(gm), 0., 0.01);
	}

	{
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);

		Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
		Gaugemomenta gm(system);
		gm.zero();

		physics::algorithms::md_update_gaugemomentum_gauge(&gm, .5, gf, system);
		BOOST_CHECK_CLOSE(squarenorm(gm), 13180.824966859615, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(md_update_gaugemomentum_fermion)
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

		pseudo_randomize<Spinorfield, spinor>(&sf1, 33);
		gm.zero();

		physics::algorithms::md_update_gaugemomentum_fermion(&gm, .6, gf, sf1, system);
		BOOST_CHECK_CLOSE(squarenorm(gm), 16956.313363328729, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(md_update_gaugemomentum_fermion_eo)
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

		pseudo_randomize<Spinorfield, spinor>(&src, 25);
		convert_to_eoprec(&sf1, &sf2, src);
		gm.zero();

		physics::algorithms::md_update_gaugemomentum_fermion(&gm, .6, gf, sf1, system);
		BOOST_CHECK_CLOSE(squarenorm(gm), 1324.8991304359643, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(md_update_gaugemomentum_detratio)
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

		pseudo_randomize<Spinorfield, spinor>(&sf1, 35);
		gm.zero();

		physics::algorithms::md_update_gaugemomentum_detratio(&gm, .4, gf, sf1, system);
		BOOST_CHECK_CLOSE(squarenorm(gm), 2.8065441263959463e-12, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(md_update_gaugemomentum_detratio_eo)
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

		pseudo_randomize<Spinorfield, spinor>(&src, 27);
		convert_to_eoprec(&sf1, &sf2, src);
		gm.zero();

		physics::algorithms::md_update_gaugemomentum_fermion(&gm, .4, gf, sf1, system);
		BOOST_CHECK_CLOSE(squarenorm(gm), 561.82056232890591, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(md_update_gaugemomentum)
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

		pseudo_randomize<Spinorfield, spinor>(&sf1, 37);
		gm.zero();

		physics::algorithms::md_update_gaugemomentum(&gm, .5, gf, sf1, system);
		BOOST_CHECK_CLOSE(squarenorm(gm), 22735.531268100596, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(md_update_gaugemomentum_eo)
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

		pseudo_randomize<Spinorfield, spinor>(&src, 29);
		convert_to_eoprec(&sf1, &sf2, src);
		gm.zero();

		physics::algorithms::md_update_gaugemomentum(&gm, .6, gf, sf1, system);
		BOOST_CHECK_CLOSE(squarenorm(gm), 20459.049258021565, 0.01);
	}
}
