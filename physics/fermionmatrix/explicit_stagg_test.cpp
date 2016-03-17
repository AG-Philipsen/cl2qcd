/** @file
 * Tests of the explicit staggered fermionmatrix implementations
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

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::fermionmatrix::explicit
#include <boost/test/unit_test.hpp>

#include "../lattices/util.hpp"
#include "../../host_functionality/logger.hpp"
#include "../test_util_staggered.h"
#include "../../hardware/code/spinors.hpp"
#include "../../interfaceImplementations/interfacesHandler.hpp"
#include "../../interfaceImplementations/hardwareParameters.hpp"
#include "../../interfaceImplementations/openClKernelParameters.hpp"

BOOST_AUTO_TEST_CASE(D_KS_eo)
{
//void physics::fermionmatrix::DKS_eo(const physics::lattices::Staggeredfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Staggeredfield_eo& in, int evenodd)
	{
		logger.info() << "First test...";
		//This test is with cold links, periodic BC, random field, 8**4 lattice
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--nspace=8", "--fermact=rooted_stagg", "--num_dev=1"};
		meta::Inputparameters params(4, _params);
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, false);
		Staggeredfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
		Staggeredfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
		Staggeredfield_eo out(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());

		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf1, 13);
		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf2, 31);
		

	//The following lines are to be used to produce the ref_vec file needed to get the ref_value
	//---> Comment them out when the reference values have been obtained!
	 /*
	print_staggeredfield_eo_to_textfile("ref_vec_odd", &sf1, system);
	print_staggeredfield_eo_to_textfile("ref_vec_even", &sf2, system);
	logger.info() << "Produced the ref_vec text file with the staggered field for the ref. code. Returning...";
	return;
	// */
		physics::fermionmatrix::DKS_eo(&out, gf, sf1, EVEN);
		BOOST_CHECK_CLOSE(squarenorm(out), 2030.1639500272767691, 1.e-8);
		physics::fermionmatrix::DKS_eo(&out, gf, sf2, ODD);
		BOOST_CHECK_CLOSE(squarenorm(out), 2076.7437224316167885, 1.e-8);
	}

	{
		logger.info() << "Second test...";
		//This test is with hot links, periodic BC, random field, 4**4 lattice
		using namespace physics::lattices;
		const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg", "--num_dev=1"};
		meta::Inputparameters params(4, _params);
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
		physics::PrngParametersImplementation prngParameters{params};
		physics::PRNG prng{system, &prngParameters};

		//This configuration for the Ref.Code is the same as for example dks_input_5
		Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
		Staggeredfield_eo sf1(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
		Staggeredfield_eo sf2(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());
		Staggeredfield_eo out(system, interfacesHandler.getInterface<physics::lattices::Staggeredfield_eo>());

		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf1, 123);
		pseudo_randomize<Staggeredfield_eo, su3vec>(&sf2, 321);

	//The following lines are to be used to produce the ref_vec file needed to get the ref_value
	//---> Comment them out when the reference values have been obtained!
	/*
	print_staggeredfield_eo_to_textfile("ref_vec_odd", &sf1, system);
	print_staggeredfield_eo_to_textfile("ref_vec_even", &sf2, system);
	logger.info() << "Produced the ref_vec text file with the staggered field for the ref. code. Returning...";
	return;
	// */

		physics::fermionmatrix::DKS_eo(&out, gf, sf1, EVEN);
		BOOST_CHECK_CLOSE(squarenorm(out), 547.69039343718509372, 1.e-8);
		physics::fermionmatrix::DKS_eo(&out, gf, sf2, ODD);
		BOOST_CHECK_CLOSE(squarenorm(out), 536.10645183266251479, 1.e-8);
	}

}

