/** @file
 * Tests of the fermion force algorithms
 *
 * Copyright (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 * Copyright (c) 2012-2013 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>
 * Copyright (c) 2013, 2017 Alessandro Sciarra <sciarra@th.phys.uni-frankfurt.de>
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

#include "forces.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::lattice::forces
#include <boost/test/unit_test.hpp>

#include "../lattices/util.hpp"
#include "../../interfaceImplementations/interfacesHandler.hpp"
#include "../../interfaceImplementations/hardwareParameters.hpp"
#include "../../interfaceImplementations/openClKernelParameters.hpp"

BOOST_AUTO_TEST_CASE(gauge_force)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=16"};
		meta::Inputparameters params(2, _params);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
        hardware::HardwareParametersImplementation hP(&params);
        hardware::code::OpenClKernelParametersImplementation kP(params);
        hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, false);
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());
		gm.zero();

		physics::algorithms::gauge_force(&gm, gf);
		BOOST_CHECK_CLOSE(squarenorm(gm), 0., 0.01);
	}

	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
        hardware::HardwareParametersImplementation hP(&params);
        hardware::code::OpenClKernelParametersImplementation kP(params);
        hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());
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
        physics::InterfacesHandlerImplementation interfacesHandler{params};
        hardware::HardwareParametersImplementation hP(&params);
        hardware::code::OpenClKernelParametersImplementation kP(params);
        hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, false);
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());
		gm.zero();

		physics::algorithms::gauge_force_tlsym(&gm, gf);
		BOOST_CHECK_CLOSE(squarenorm(gm), 0., 0.01);
	}

	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4", "--gaugeact=tlsym"};
		meta::Inputparameters params(3, _params);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
        hardware::HardwareParametersImplementation hP(&params);
        hardware::code::OpenClKernelParametersImplementation kP(params);
        hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());
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
        physics::InterfacesHandlerImplementation interfacesHandler{params};
        hardware::HardwareParametersImplementation hP(&params);
        hardware::code::OpenClKernelParametersImplementation kP(params);
        hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, false);
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());
		gm.zero();

		physics::algorithms::calc_gauge_force(&gm, gf, interfacesHandler);
		BOOST_CHECK_CLOSE(squarenorm(gm), 0., 0.01);
	}

	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
        hardware::HardwareParametersImplementation hP(&params);
        hardware::code::OpenClKernelParametersImplementation kP(params);
        hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());
		gm.zero();

		physics::algorithms::calc_gauge_force(&gm, gf, interfacesHandler);
		BOOST_CHECK_CLOSE(squarenorm(gm), 52723.299867438458, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(calc_tot_force)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
        hardware::HardwareParametersImplementation hP(&params);
        hardware::code::OpenClKernelParametersImplementation kP(params);
        hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Spinorfield sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());

		pseudo_randomize<Spinorfield, spinor>(&sf1, 24);
		gm.zero();

		physics::algorithms::calc_total_force(&gm, gf, sf1, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Spinorfield>());
		BOOST_CHECK_CLOSE(squarenorm(gm), 89317.106966900712, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(calc_tot_force_eo)
{
	{
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4"};
		meta::Inputparameters params(2, _params);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
        hardware::HardwareParametersImplementation hP(&params);
        hardware::code::OpenClKernelParametersImplementation kP(params);
        hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Spinorfield src(system, interfacesHandler.getInterface<physics::lattices::Spinorfield>());
		Spinorfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
		Spinorfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Spinorfield_eo>());
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());

		pseudo_randomize<Spinorfield, spinor>(&src, 25);
		convert_to_eoprec(&sf1, &sf2, src);
		gm.zero();

		physics::algorithms::calc_total_force(&gm, gf, sf1, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Spinorfield_eo>());
		BOOST_CHECK_CLOSE(squarenorm(gm), 56762.555327447422, 0.01);
	}
}

BOOST_AUTO_TEST_CASE(calc_tot_stagg_force_eo)
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
		Rooted_Staggeredfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Rooted_Staggeredfield_eo>(), approx);
		Rooted_Staggeredfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Rooted_Staggeredfield_eo>(), approx);
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());
		
		//These are the same fields of the excplicit test D_KS_eo (second test)
		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf1[0], 123); //it will be A
		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf2[0], 321); //it will be B
		
		gm.zero();
		calc_total_force(&gm, gf, sf1, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Rooted_Staggeredfield_eo>());
		BOOST_CHECK_CLOSE(squarenorm(gm), 58639.680325374203676, 1.e-6);
		
		gm.zero();
		calc_total_force(&gm, gf, sf2, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Rooted_Staggeredfield_eo>());
		BOOST_CHECK_CLOSE(squarenorm(gm), 57864.102469501536689, 1.e-6);
	}
	
	{
		const char * _params[] = {"foo", "--ntime=4", "--theta_fermion_temporal=1", "--fermact=rooted_stagg", "--mass=0.125", "--num_dev=1"};
		meta::Inputparameters params(6, _params);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
        hardware::HardwareParametersImplementation hP(&params);
        hardware::code::OpenClKernelParametersImplementation kP(params);
        hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};
		
		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Rooted_Staggeredfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Rooted_Staggeredfield_eo>(), approx);
		Rooted_Staggeredfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Rooted_Staggeredfield_eo>(), approx);
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());
		
		//These are the same fields of the excplicit test D_KS_eo (second test)
		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf1[0], 123); //it will be A
		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf2[0], 321); //it will be B
		
		gm.zero();
		calc_total_force(&gm, gf, sf1, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Rooted_Staggeredfield_eo>());
		BOOST_CHECK_CLOSE(squarenorm(gm), 58424.656915726853185, 1.e-6);
		
		gm.zero();
		calc_total_force(&gm, gf, sf2, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Rooted_Staggeredfield_eo>());
		BOOST_CHECK_CLOSE(squarenorm(gm), 58492.589653369606822, 1.e-6);
	}
	
	{
		const char * _params[] = {"foo", "--ntime=4", "--theta_fermion_temporal=1", "--fermact=rooted_stagg",
		                          "--use_chem_pot_im=true", "--chem_pot_im=0.5678", "--mass=0.125", "--num_dev=1"};
		meta::Inputparameters params(8, _params);
        physics::InterfacesHandlerImplementation interfacesHandler{params};
        hardware::HardwareParametersImplementation hP(&params);
        hardware::code::OpenClKernelParametersImplementation kP(params);
        hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};
		
		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Rooted_Staggeredfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Rooted_Staggeredfield_eo>(), approx);
		Rooted_Staggeredfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Rooted_Staggeredfield_eo>(), approx);
		Gaugemomenta gm(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());
		
		//These are the same fields of the excplicit test D_KS_eo (second test)
		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf1[0], 123); //it will be A
		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf2[0], 321); //it will be B
		
		gm.zero();
		calc_total_force(&gm, gf, sf1, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Rooted_Staggeredfield_eo>());
		BOOST_CHECK_CLOSE(squarenorm(gm), 57415.997451495910354, 1.e-6);
		
		gm.zero();
		calc_total_force(&gm, gf, sf2, system, interfacesHandler, interfacesHandler.getAdditionalParameters<Rooted_Staggeredfield_eo>());
		BOOST_CHECK_CLOSE(squarenorm(gm), 57338.140878283666098, 1.e-6);
	}
	
}
