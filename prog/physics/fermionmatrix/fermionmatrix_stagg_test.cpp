/** @file
 * Tests of the fermionmatrix implementations
 * 
 * (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
 */

#include "fermionmatrix_stagg.hpp"

#include "../lattices/util.hpp"
#include <boost/type_traits.hpp>
#include <boost/utility.hpp>

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::fermionmatrix::<fermionmatrix_stagg>
#include <boost/test/unit_test.hpp>

//So far, only a test for even-odd fermionmatrix object is needed
/*
template<class FERMIONMATRIX>
typename boost::enable_if<boost::is_base_of<physics::fermionmatrix::Fermionmatrix, FERMIONMATRIX>, void>::type
test_fermionmatrix_stagg(const hmc_float refs[4], const int seed);*/
template<class FERMIONMATRIX> typename boost::enable_if<boost::is_base_of<physics::fermionmatrix::Fermionmatrix_stagg_eo, FERMIONMATRIX>, void>::type test_fermionmatrix_stagg(const hmc_float refs[4]);
//Specialization of the template above
template<> typename boost::enable_if<boost::is_base_of<physics::fermionmatrix::Fermionmatrix_stagg_eo, physics::fermionmatrix::D_KS_eo>, void>::type test_fermionmatrix_stagg <physics::fermionmatrix::D_KS_eo>(const hmc_float refs[4]);


BOOST_AUTO_TEST_CASE(D_KS_eo)
{
	const hmc_float refs[4] = {2030.1639500272767691, 2076.7437224316167885,
	                           547.69039343718509372, 536.10645183266251479};
	test_fermionmatrix_stagg<physics::fermionmatrix::D_KS_eo>(refs);
}

BOOST_AUTO_TEST_CASE(MdagM_eo)
{
	const hmc_float refs[4] = {12670.438549218517437, 12925.274439556447760,
	                           2980.0540282922056576, 2883.0056841550949684};
	test_fermionmatrix_stagg<physics::fermionmatrix::MdagM_eo>(refs);
}


template<class FERMIONMATRIX>
typename boost::enable_if<boost::is_base_of<physics::fermionmatrix::Fermionmatrix_stagg_eo, FERMIONMATRIX>, void>::type test_fermionmatrix_stagg(const hmc_float refs[4])
{
	{
		//Test with cold links, periodic BC, random field, 8**4 lattice
		logger.info() << "First test...";
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--nspace=8"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);
		
		FERMIONMATRIX matrix1(system, 1., ODD);
		FERMIONMATRIX matrix2(system, 1.);

		logger.info() << "The mass of the fermion is " << matrix1.get_mass();
		
		Gaugefield gf(system, prng, false);
		Staggeredfield_eo sf1(system);
		Staggeredfield_eo sf2(system);
		Staggeredfield_eo out(system);

		//Using these seeds I can use for the reference code the fermionic 
		//fields of the explicit_stagg_test file, but pay attention to the fact
		//that in explicit_stagg_test sf1 is an odd field: you must be coherent!
		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf1, 13);
		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf2, 31);

		matrix1(&out, gf, sf1);
		BOOST_CHECK_CLOSE(squarenorm(out), refs[0], 1.e-8);
		matrix2(&out, gf, sf2);
		BOOST_CHECK_CLOSE(squarenorm(out), refs[1], 1.e-8);
	}

	{
		//Test with hot links, periodic BC, random field, 4**4 lattice
		logger.info() << "Second test...";
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);
		
		FERMIONMATRIX matrix1(system, 1., ODD);
		FERMIONMATRIX matrix2(system, 1.);

		//This configuration for the Ref.Code is the same as for example dks_input_5
		Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
		Staggeredfield_eo sf1(system);
		Staggeredfield_eo sf2(system);
		Staggeredfield_eo out(system);

		//Using these seeds I can use for the reference code the fermionic 
		//fields of the explicit_stagg_test file, but pay attention to the fact
		//that in explicit_stagg_test sf1 is an odd field: you must be coherent!
		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf1, 123);
		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf2, 321);

		matrix1(&out, gf, sf1);
		BOOST_CHECK_CLOSE(squarenorm(out), refs[2], 1.e-8);
		matrix2(&out, gf, sf2);
		BOOST_CHECK_CLOSE(squarenorm(out), refs[3], 1.e-8);
	}
}


template<> typename boost::enable_if<boost::is_base_of<physics::fermionmatrix::Fermionmatrix_stagg_eo, physics::fermionmatrix::D_KS_eo>, void>::type test_fermionmatrix_stagg <physics::fermionmatrix::D_KS_eo>(const hmc_float refs[4])
{
	{
		//Test with cold links, periodic BC, random field, 8**4 lattice
		logger.info() << "First test...";
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--nspace=8"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);
		
		physics::fermionmatrix::D_KS_eo matrix1(system, EVEN);
		physics::fermionmatrix::D_KS_eo matrix2(system, ODD);

		Gaugefield gf(system, prng, false);
		Staggeredfield_eo sf1(system);
		Staggeredfield_eo sf2(system);
		Staggeredfield_eo out(system);

		//Using these seeds I can use for the reference code the fermionic 
		//fields of the explicit_stagg_test file, but pay attention to the fact
		//that in explicit_stagg_test sf1 is an odd field: you must be coherent!
		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf1, 13);
		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf2, 31);

		matrix1(&out, gf, sf1);
		BOOST_CHECK_CLOSE(squarenorm(out), refs[0], 1.e-8);
		matrix2(&out, gf, sf2);
		BOOST_CHECK_CLOSE(squarenorm(out), refs[1], 1.e-8);
	}

	{
		//Test with hot links, periodic BC, random field, 4**4 lattice
		logger.info() << "Second test...";
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);
		
		physics::fermionmatrix::D_KS_eo matrix1(system, EVEN);
		physics::fermionmatrix::D_KS_eo matrix2(system, ODD);

		//This configuration for the Ref.Code is the same as for example dks_input_5
		Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
		Staggeredfield_eo sf1(system);
		Staggeredfield_eo sf2(system);
		Staggeredfield_eo out(system);

		//Using these seeds I can use for the reference code the fermionic 
		//fields of the explicit_stagg_test file, but pay attention to the fact
		//that in explicit_stagg_test sf1 is an odd field: you must be coherent!
		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf1, 123);
		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf2, 321);

		matrix1(&out, gf, sf1);
		BOOST_CHECK_CLOSE(squarenorm(out), refs[2], 1.e-8);
		matrix2(&out, gf, sf2);
		BOOST_CHECK_CLOSE(squarenorm(out), refs[3], 1.e-8);
	}
}




/*
template<class FERMIONMATRIX>
typename boost::enable_if<boost::is_base_of<physics::fermionmatrix::Fermionmatrix, FERMIONMATRIX>, void>::type
test_fermionmatrix(const hmc_float refs[4], const int seed)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=16"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);
		FERMIONMATRIX matrix(ARG_DEF, ARG_DEF, system);

		Gaugefield gf(system, prng, false);
		Spinorfield sf1(system);
		Spinorfield sf2(system);

		pseudo_randomize<Spinorfield, spinor>(&sf1, seed);
		pseudo_randomize<Spinorfield, spinor>(&sf2, seed + 1);

		matrix(&sf2, gf, sf1);
		BOOST_CHECK_CLOSE(squarenorm(sf2), refs[0], 0.01);
		matrix(&sf1, gf, sf2);
		BOOST_CHECK_CLOSE(squarenorm(sf1), refs[1], 0.01);
	}

	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
		hardware::System system(params);
		physics::PRNG prng(system);
		FERMIONMATRIX matrix(ARG_DEF, ARG_DEF, system);

		Gaugefield gf(system, prng, std::string(SOURCEDIR) + "/tests/conf.00200");
		Spinorfield sf1(system);
		Spinorfield sf2(system);

		pseudo_randomize<Spinorfield, spinor>(&sf1, seed + 2);
		pseudo_randomize<Spinorfield, spinor>(&sf2, seed + 3);

		matrix(&sf2, gf, sf1);
		BOOST_CHECK_CLOSE(squarenorm(sf2), refs[2], 0.01);
		matrix(&sf1, gf, sf2);
		BOOST_CHECK_CLOSE(squarenorm(sf1), refs[3], 0.01);
	}
}
*/

