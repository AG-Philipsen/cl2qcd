/** @file
 * Tests of the fermion force algorithms
 */

#include "forces.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::forces
#include <boost/test/unit_test.hpp>

#include "../lattices/util.hpp"

BOOST_AUTO_TEST_CASE(gauge_force)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=16"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);

		Gaugefield gf(system, prng, false);
		Gaugemomenta gm(system);
		gm.zero();

		physics::algorithms::gauge_force(&gm, gf);
		BOOST_CHECK_CLOSE(squarenorm(gm), 0., 0.01);
	}

	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);

		Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
		Gaugemomenta gm(system);
		gm.zero();

		physics::algorithms::gauge_force(&gm, gf);
		BOOST_CHECK_CLOSE(squarenorm(gm), 52723.299867438458, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(gauge_force_tlsym)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=16", "--gaugeact=tlsym"};
		meta::Inputparameters params(3, _params);
		hardware::System system(params);
		physics::PRNG prng(system);

		Gaugefield gf(system, prng, false);
		Gaugemomenta gm(system);
		gm.zero();

		physics::algorithms::gauge_force_tlsym(&gm, gf);
		BOOST_CHECK_CLOSE(squarenorm(gm), 0., 0.01);
	}

	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4", "--gaugeact=tlsym"};
		meta::Inputparameters params(3, _params);
		hardware::System system(params);
		physics::PRNG prng(system);

		Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
		Gaugemomenta gm(system);
		gm.zero();

		physics::algorithms::gauge_force_tlsym(&gm, gf);
		BOOST_CHECK_CLOSE(squarenorm(gm), 2016.6355154119337, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(calc_gauge_force)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=16"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);

		Gaugefield gf(system, prng, false);
		Gaugemomenta gm(system);
		gm.zero();

		physics::algorithms::calc_gauge_force(&gm, gf, system);
		BOOST_CHECK_CLOSE(squarenorm(gm), 0., 0.01);
	}

	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);

		Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
		Gaugemomenta gm(system);
		gm.zero();

		physics::algorithms::calc_gauge_force(&gm, gf, system);
		BOOST_CHECK_CLOSE(squarenorm(gm), 52723.299867438458, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(calc_tot_force)
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

		pseudo_randomize<Spinorfield, spinor>(&sf1, 24);
		gm.zero();

		physics::algorithms::calc_total_force(&gm, gf, sf1, system);
		BOOST_CHECK_CLOSE(squarenorm(gm), 89317.106966900712, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(calc_tot_force_eo)
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

		physics::algorithms::calc_total_force(&gm, gf, sf1, system);
		BOOST_CHECK_CLOSE(squarenorm(gm), 56762.555327447422, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(calc_tot_stagg_force_eo)
{
	using namespace physics::lattices;
	using namespace physics::algorithms;
	Rational_Approximation approx(8, 1,2, 1.e-5,1);
	{
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);
		
		Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
		Rooted_Staggeredfield_eo sf1(approx, system);
		Rooted_Staggeredfield_eo sf2(approx, system);
		Gaugemomenta gm(system);
		
		//These are the same fields of the excplicit test D_KS_eo (second test)
		pseudo_randomize<Rooted_Staggeredfield_eo, su3vec>(&sf1, 123); //it will be A
		pseudo_randomize<Rooted_Staggeredfield_eo, su3vec>(&sf2, 321); //it will be B
		
		gm.zero();
		calc_total_force(&gm, gf, sf1, system, 0.125);
		BOOST_CHECK_CLOSE(squarenorm(gm), 58639.680325374203676, 1.e-6);
		
		gm.zero();
		calc_total_force(&gm, gf, sf2, system, 0.125);
		BOOST_CHECK_CLOSE(squarenorm(gm), 57864.102469501536689, 1.e-6);
	}
	
	{
		const char * _params[] = {"foo", "--ntime=4", "--theta_fermion_temporal=1"};
		meta::Inputparameters params(3, _params);
		hardware::System system(params);
		physics::PRNG prng(system);
		
		Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
		Rooted_Staggeredfield_eo sf1(approx, system);
		Rooted_Staggeredfield_eo sf2(approx, system);
		Gaugemomenta gm(system);
		
		//These are the same fields of the excplicit test D_KS_eo (second test)
		pseudo_randomize<Rooted_Staggeredfield_eo, su3vec>(&sf1, 123); //it will be A
		pseudo_randomize<Rooted_Staggeredfield_eo, su3vec>(&sf2, 321); //it will be B
		
		gm.zero();
		calc_total_force(&gm, gf, sf1, system, 0.125);
		BOOST_CHECK_CLOSE(squarenorm(gm), 58424.656915726853185, 1.e-6);
		
		gm.zero();
		calc_total_force(&gm, gf, sf2, system, 0.125);
		BOOST_CHECK_CLOSE(squarenorm(gm), 58492.589653369606822, 1.e-6);
	}
	
}
