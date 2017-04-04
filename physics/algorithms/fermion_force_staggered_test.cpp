/** @file
 * Tests of the fermion force algorithms
 *
 * Copyright (c) 2013, 2017 Alessandro Sciarra <sciarra@th.phys.uni-frankfurt.de>
 * Copyright (c) 2017 Francesca Cuteri <cuteri@th.physik.uni-frankfurt.de>
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

#include "fermion_force_staggered.hpp"
#include "../../interfaceImplementations/interfacesHandler.hpp"
#include "../../interfaceImplementations/hardwareParameters.hpp"
#include "../../interfaceImplementations/openClKernelParameters.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::fermion_force
#include <boost/test/unit_test.hpp>


BOOST_AUTO_TEST_CASE(fermion_force_staggered_eo)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg", "--num_dev=1"};
		meta::Inputparameters params(4, _params);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
        hardware::HardwareParametersImplementation hP(&params);
        hardware::code::OpenClKernelParametersImplementation kP(params);
        hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, false);
		Staggeredfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
		Staggeredfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());
		
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
		const char * _params[] = {"foo", "--ntime=16", "--nspace=8", "--fermact=rooted_stagg", "--num_dev=1"};
		meta::Inputparameters params(5, _params);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
        hardware::HardwareParametersImplementation hP(&params);
        hardware::code::OpenClKernelParametersImplementation kP(params);
        hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};
		
		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, false);
		Staggeredfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
		Staggeredfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());

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
		const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg", "--num_dev=1"};
		meta::Inputparameters params(4, _params);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
        hardware::HardwareParametersImplementation hP(&params);
        hardware::code::OpenClKernelParametersImplementation kP(params);
        hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Staggeredfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
		Staggeredfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());
		
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
	
	//These are with imaginary chemical potential different from zero
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--nspace=8", "--fermact=rooted_stagg", "--use_chem_pot_im=true", "--chem_pot_im=0.5678", "--num_dev=1"};
		meta::Inputparameters params(6, _params);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
        hardware::HardwareParametersImplementation hP(&params);
        hardware::code::OpenClKernelParametersImplementation kP(params);
        hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};
		
		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, false);
		Staggeredfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
		Staggeredfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());

		sf1.set_cold();
		sf2.set_cold();
		gm.zero();
		
		//If links are one and fermionic fields are set to cold, one can analitically
		//calculate the result: it is 2*sin(chem_pot_im)^2 / (3 * vol^2) for each time direction
		//for the norm of the su3 matrix. Then here gm is the 8 component vector with the
		//algebra coefficients and then we have still a factor 2 in the squarenorm.
		BOOST_REQUIRE_EQUAL(squarenorm(gm), 0);
		physics::algorithms::fermion_force(&gm, sf1, sf2, gf, EVEN);
		//Here, for example: 2*[2*sin(chem_pot_im)^2 / (3 * vol^2)] * vol/2
		BOOST_REQUIRE_CLOSE(squarenorm(gm), 0.0000470712535699, 1.e-8);
		physics::algorithms::fermion_force(&gm, sf1, sf2, gf, ODD);
		BOOST_REQUIRE_CLOSE(squarenorm(gm), 0.0000941425071398, 1.e-8);
	}
	
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg", "--use_chem_pot_im=true", "--chem_pot_im=0.5678", "--theta_fermion_temporal=1.", "--num_dev=1"};
		meta::Inputparameters params(7, _params);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
        hardware::HardwareParametersImplementation hP(&params);
        hardware::code::OpenClKernelParametersImplementation kP(params);
        hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Staggeredfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
		Staggeredfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());
		
		//These are the same fields of the excplicit test D_KS_eo (second test)
		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf1, 123); //it will be A
		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf2, 321); //it will be B
		gm.zero();

		BOOST_REQUIRE_EQUAL(squarenorm(gm), 0);
		physics::algorithms::fermion_force(&gm, sf1, sf2, gf, EVEN);
		BOOST_CHECK_CLOSE(squarenorm(gm), 1950.7859609652412018, 1.e-8);
		//Note that now the ODD part is added to the EVEN one
		physics::algorithms::fermion_force(&gm, sf1, sf2, gf, ODD);
		BOOST_CHECK_CLOSE(squarenorm(gm), 3919.7908258193665461, 1.e-8);
	}
	
}


BOOST_AUTO_TEST_CASE(calc_fermion_force_staggered_eo)
{
	using namespace physics::lattices;
	using namespace physics::algorithms;
	
	Rational_Approximation approx(8, 1,2, 1.e-5,1);
	
	{
		const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg", "--mass=0.125", "--num_dev=1"};
		meta::Inputparameters params(5, _params);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());
		Rooted_Staggeredfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Rooted_Staggeredfield_eo>(), approx);
		Rooted_Staggeredfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Rooted_Staggeredfield_eo>(), approx);
	
		//These are the same fields of the excplicit test D_KS_eo (second test)
		pseudo_randomize<Staggeredfield_eo, su3vec>(sf1[0].get(), 123); //it will be A
		pseudo_randomize<Staggeredfield_eo, su3vec>(sf2[0].get(), 321); //it will be B
	
		gm.zero();
		physics::algorithms::calc_fermion_forces(&gm, gf, sf1, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Rooted_Staggeredfield_eo>());
		BOOST_CHECK_CLOSE(squarenorm(gm), 2214.9003939576623452, 1.e-6);
	
		gm.zero();
		physics::algorithms::calc_fermion_forces(&gm, gf, sf2, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Rooted_Staggeredfield_eo>());
		BOOST_CHECK_CLOSE(squarenorm(gm), 1845.6513833002247793, 1.e-6);
	}
	
	{
		using namespace physics::lattices;
		using namespace physics::algorithms;
		const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg", "--use_chem_pot_im=true", "--mass=0.125", "--chem_pot_im=0.5678", "--num_dev=1"};
		meta::Inputparameters params(7, _params);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
        hardware::HardwareParametersImplementation hP(&params);
        hardware::code::OpenClKernelParametersImplementation kP(params);
        hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());
		Rooted_Staggeredfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Rooted_Staggeredfield_eo>(), approx);
		Rooted_Staggeredfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Rooted_Staggeredfield_eo>(), approx);
	
		//These are the same fields of the excplicit test D_KS_eo (second test)
		pseudo_randomize<Staggeredfield_eo, su3vec>(sf1[0].get(), 123); //it will be A
		pseudo_randomize<Staggeredfield_eo, su3vec>(sf2[0].get(), 321); //it will be B
	
		gm.zero();
		physics::algorithms::calc_fermion_forces(&gm, gf, sf1, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Rooted_Staggeredfield_eo>());
		BOOST_CHECK_CLOSE(squarenorm(gm), 2315.4175592593164765, 1.e-6);
	
		gm.zero();
		physics::algorithms::calc_fermion_forces(&gm, gf, sf2, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Rooted_Staggeredfield_eo>());
		BOOST_CHECK_CLOSE(squarenorm(gm), 1932.6440761489629949, 1.e-6);
	}	
}












