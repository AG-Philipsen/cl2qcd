/** @file
 * Tests of the fermionmatrix implementations
 *
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "fermionmatrix_stagg.hpp"

#include "../../interfaceImplementations/interfacesHandler.hpp"
#include "../../interfaceImplementations/hardwareParameters.hpp"
#include "../../interfaceImplementations/openClKernelParameters.hpp"
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
		const char * _params[] = {"foo", "--nspace=8", "--fermact=rooted_stagg", "--mass=1.", "--num_dev=1"};
		meta::Inputparameters params(5, _params);
		GaugefieldParametersImplementation gaugefieldParameters( &params );
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};
		
		FERMIONMATRIX matrix1(system, interfacesHandler.getInterface<FERMIONMATRIX>(), ODD);
		FERMIONMATRIX matrix2(system, interfacesHandler.getInterface<FERMIONMATRIX>());

		logger.info() << "The mass of the fermion is " << params.get_mass();
		
		Gaugefield gf(system, &gaugefieldParameters, prng, false);
		Staggeredfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
		Staggeredfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
		Staggeredfield_eo out(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());

		//Using these seeds I can use for the reference code the fermionic 
		//fields of the explicit_stagg_test file, but pay attention to the fact
		//that in explicit_stagg_test sf1 is an odd field: you must be coherent!
		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf1, 13);
		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf2, 31);

		matrix1(&out, gf, sf1, &(interfacesHandler.getAdditionalParameters<Staggeredfield_eo>()));
		BOOST_CHECK_CLOSE(squarenorm(out), refs[0], 1.e-8);
		matrix2(&out, gf, sf2, &(interfacesHandler.getAdditionalParameters<Staggeredfield_eo>()));
		BOOST_CHECK_CLOSE(squarenorm(out), refs[1], 1.e-8);
	}

	{
		//Test with hot links, periodic BC, random field, 4**4 lattice
		logger.info() << "Second test...";
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg", "--mass=1.", "--num_dev=1"};
		meta::Inputparameters params(5, _params);
		GaugefieldParametersImplementation gaugefieldParameters( &params );
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};
		
		FERMIONMATRIX matrix1(system, interfacesHandler.getInterface<FERMIONMATRIX>(), ODD);
		FERMIONMATRIX matrix2(system, interfacesHandler.getInterface<FERMIONMATRIX>());

		//This configuration for the Ref.Code is the same as for example dks_input_5
		Gaugefield gf(system, &gaugefieldParameters, prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Staggeredfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
		Staggeredfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
		Staggeredfield_eo out(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());

		//Using these seeds I can use for the reference code the fermionic 
		//fields of the explicit_stagg_test file, but pay attention to the fact
		//that in explicit_stagg_test sf1 is an odd field: you must be coherent!
		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf1, 123);
		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf2, 321);

		matrix1(&out, gf, sf1, &(interfacesHandler.getAdditionalParameters<Staggeredfield_eo>()));
		BOOST_CHECK_CLOSE(squarenorm(out), refs[2], 1.e-8);
		matrix2(&out, gf, sf2, &(interfacesHandler.getAdditionalParameters<Staggeredfield_eo>()));
		BOOST_CHECK_CLOSE(squarenorm(out), refs[3], 1.e-8);
	}
}


template<> typename boost::enable_if<boost::is_base_of<physics::fermionmatrix::Fermionmatrix_stagg_eo, physics::fermionmatrix::D_KS_eo>, void>::type test_fermionmatrix_stagg <physics::fermionmatrix::D_KS_eo>(const hmc_float refs[4])
{
	{
		//Test with cold links, periodic BC, random field, 8**4 lattice
		logger.info() << "First test...";
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--nspace=8", "--fermact=rooted_stagg", "--num_dev=1"};
		meta::Inputparameters params(4, _params);
		GaugefieldParametersImplementation gaugefieldParameters( &params );
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};
		
		physics::fermionmatrix::D_KS_eo matrix1(system, interfacesHandler.getInterface<physics::fermionmatrix::D_KS_eo>(), EVEN);
		physics::fermionmatrix::D_KS_eo matrix2(system, interfacesHandler.getInterface<physics::fermionmatrix::D_KS_eo>(), ODD);

		Gaugefield gf(system, &gaugefieldParameters, prng, false);
		Staggeredfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
		Staggeredfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
		Staggeredfield_eo out(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());

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
		const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg", "--num_dev=1"};
		meta::Inputparameters params(4, _params);
		GaugefieldParametersImplementation gaugefieldParameters( &params );
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};
		
		physics::fermionmatrix::D_KS_eo matrix1(system, interfacesHandler.getInterface<physics::fermionmatrix::D_KS_eo>(), EVEN);
		physics::fermionmatrix::D_KS_eo matrix2(system, interfacesHandler.getInterface<physics::fermionmatrix::D_KS_eo>(), ODD);

		//This configuration for the Ref.Code is the same as for example dks_input_5
		Gaugefield gf(system, &gaugefieldParameters, prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Staggeredfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
		Staggeredfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
		Staggeredfield_eo out(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());

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


