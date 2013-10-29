/** @file
 * Tests of the fermion force algorithms
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 * (c) 2013 Alessandro Sciarra <sciarra@compeng.uni-frankfurt.de>
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

		BOOST_REQUIRE_SMALL(squarenorm(gm), 0.001);
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
		const char * _params[] = {"foo", "--ntime=4", "--kappa_mp=.25"};
		meta::Inputparameters params(3, _params);
		hardware::System system(params);
		physics::PRNG prng(system);

		Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
		Spinorfield sf1(system);
		Gaugemomenta gm(system);

		pseudo_randomize<Spinorfield, spinor>(&sf1, 21);
		gm.zero();

		physics::algorithms::calc_detratio_forces(&gm, gf, sf1, system);
		BOOST_CHECK_CLOSE(squarenorm(gm), 12480.139807156647, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(calc_detratio_forces_eo)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4", "--kappa_mp=.25"};
		meta::Inputparameters params(3, _params);
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
		BOOST_CHECK_CLOSE(squarenorm(gm), 33313.511647643441, 0.01);
	}
}

////////////////////////////////////////////////////////////////////////////
// STAGGERED TESTS

BOOST_AUTO_TEST_CASE(fermion_force_staggered_eo)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);

		Gaugefield gf(system, prng, false);
		Staggeredfield_eo sf1(system);
		Staggeredfield_eo sf2(system);
		Gaugemomenta gm(system);
		
		//These are the same fields of the excplicit test D_KS_eo (second test)
		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf1, 123); //it will be A
		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf2, 321); //it will be B
		gm.zero();

		BOOST_REQUIRE_EQUAL(squarenorm(gm), 0);
		physics::algorithms::fermion_force(&gm, sf1, sf2, gf, EVEN);
		BOOST_CHECK_CLOSE(squarenorm(gm), 855.08060572822057566, 1.e-8);
		//Note that now the ODD part is added to the EVEN one
		physics::algorithms::fermion_force(&gm, sf1, sf2, gf, ODD);
		BOOST_CHECK_CLOSE(squarenorm(gm), 1714.8417937241449636, 1.e-8);
	}
	
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=16", "--nspace=8"};
		meta::Inputparameters params(3, _params);
		hardware::System system(params);
		physics::PRNG prng(system);
		
		Gaugefield gf(system, prng, false);
		Staggeredfield_eo sf1(system);
		Staggeredfield_eo sf2(system);
		Gaugemomenta gm(system);

		sf1.set_cold();
		sf2.set_cold();
		gm.zero();

		BOOST_REQUIRE_EQUAL(squarenorm(gm), 0);
		physics::algorithms::fermion_force(&gm, sf1, sf2, gf, EVEN);
		BOOST_REQUIRE_EQUAL(squarenorm(gm), 0);
		physics::algorithms::fermion_force(&gm, sf1, sf2, gf, ODD);
		BOOST_REQUIRE_EQUAL(squarenorm(gm), 0);
	}
	
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);

		Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
		Staggeredfield_eo sf1(system);
		Staggeredfield_eo sf2(system);
		Gaugemomenta gm(system);
		
		//These are the same fields of the excplicit test D_KS_eo (second test)
		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf1, 123); //it will be A
		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf2, 321); //it will be B
		gm.zero();

		BOOST_REQUIRE_EQUAL(squarenorm(gm), 0);
		physics::algorithms::fermion_force(&gm, sf1, sf2, gf, EVEN);
		BOOST_CHECK_CLOSE(squarenorm(gm), 1995.1105623150542669, 1.e-8);
		//Note that now the ODD part is added to the EVEN one
		physics::algorithms::fermion_force(&gm, sf1, sf2, gf, ODD);
		BOOST_CHECK_CLOSE(squarenorm(gm), 3977.231580060397846, 1.e-8);
	}
}


BOOST_AUTO_TEST_CASE(calc_fermion_force_staggered_eo)
{
	using namespace physics::lattices;
	using namespace physics::algorithms;
	const char * _params[] = {"foo", "--ntime=4"};
	meta::Inputparameters params(2, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
	Gaugemomenta gm(system);
	Rational_Approximation approx(8, 1,2, 1.e-5,1);
	Rooted_Staggeredfield_eo sf1(approx, system);
	Rooted_Staggeredfield_eo sf2(approx, system);
	
	//These are the same fields of the excplicit test D_KS_eo (second test)
	pseudo_randomize<Staggeredfield_eo, su3vec>(&sf1, 123); //it will be A
	pseudo_randomize<Staggeredfield_eo, su3vec>(&sf2, 321); //it will be B
	
	gm.zero();
	physics::algorithms::calc_fermion_forces(&gm, gf, sf1, system, 0.125);
	BOOST_CHECK_CLOSE(squarenorm(gm), 2214.9003939576623452, 1.e-6);
	
	gm.zero();
	physics::algorithms::calc_fermion_forces(&gm, gf, sf2, system, 0.125);
	BOOST_CHECK_CLOSE(squarenorm(gm), 1845.6513833002247793, 1.e-6);
}












