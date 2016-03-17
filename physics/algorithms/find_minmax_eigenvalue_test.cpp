/** @file
 * Implementation of the algorithm to find min and max eigenvalues of an operator
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

#include "find_minmax_eigenvalue.hpp"
#include "../lattices/staggeredfield_eo.hpp"
#include "../lattices/util.hpp"
#include "../../host_functionality/logger.hpp"
#include "../../interfaceImplementations/interfacesHandler.hpp"
#include "../../interfaceImplementations/hardwareParameters.hpp"
#include "../../interfaceImplementations/openClKernelParameters.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::algorithms::find_minmax_eigenvalue
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(max)
{
	using namespace physics::lattices;
	using namespace physics::algorithms;
	
	hmc_float ref_max_eig = 5.9887256245527069609;
	
	const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg", "--mass=1.01335", "--num_dev=1"};
	meta::Inputparameters params(5, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	physics::PrngParametersImplementation prngParameters{params};
	physics::PRNG prng{system, &prngParameters};
	
	//Operator for the test
	physics::fermionmatrix::MdagM_eo matrix(system, interfacesHandler.getInterface<physics::fermionmatrix::MdagM_eo>());
	//This configuration for the Ref.Code is the same as for example dks_input_5
	Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
	
	hmc_float max = find_max_eigenvalue(matrix, gf, system, interfacesHandler, 1.e-5, interfacesHandler.getAdditionalParameters<Staggeredfield_eo>());
	
	logger.info() << "ref_max_eig = " << std::setprecision(16) << ref_max_eig;
	logger.info() << "    max_eig = " << std::setprecision(16) << max;
	
	//The precision of this test is not 1.e-8 because the method is not exactly
	//in the same way implemented in the Ref.Code
	BOOST_CHECK_CLOSE(ref_max_eig, max, 1.e-6);
}


BOOST_AUTO_TEST_CASE(min)
{
	using namespace physics::lattices;
	using namespace physics::algorithms;
	
	hmc_float ref_min_eig = 1.053927941164244;
	
	const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg", "--mass=1.01335", "--num_dev=1"};
	meta::Inputparameters params(5, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	physics::PrngParametersImplementation prngParameters{params};
	physics::PRNG prng{system, &prngParameters};
	
	//Operator for the test
	physics::fermionmatrix::MdagM_eo matrix(system, interfacesHandler.getInterface<physics::fermionmatrix::MdagM_eo>());
	//This configuration for the Ref.Code is the same as for example dks_input_5
	Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
	
	hmc_float min = find_min_eigenvalue(matrix, gf, system, interfacesHandler, 1.e-3, interfacesHandler.getAdditionalParameters<Staggeredfield_eo>());
	
	logger.info() << "mass squared = " << std::setprecision(16) << params.get_mass()*params.get_mass();
	logger.info() << " ref_min_eig = " << std::setprecision(16) << ref_min_eig;
	logger.info() << "     min_eig = " << std::setprecision(16) << min;
	
	BOOST_REQUIRE_SMALL(params.get_mass()*params.get_mass(), min);
	//The precision of this test is not so high because the method to calculate
	// the eigenvalue is not exactly the same as that implemented in the Ref.Code
	BOOST_CHECK_CLOSE(ref_min_eig, min, 1.e-3);
}


BOOST_AUTO_TEST_CASE(maxmin)
{
	using namespace physics::lattices;
	using namespace physics::algorithms;
	
	hmc_float ref_max_eig = 5.2827838704124030;
	hmc_float ref_min_eig = 0.3485295092571166;
	
	const char * _params[] = {"foo", "--ntime=4", "--fermact=rooted_stagg", "--mass=0.567", "--num_dev=1"};
	meta::Inputparameters params(5, _params);
    hardware::HardwareParametersImplementation hP(&params);
    hardware::code::OpenClKernelParametersImplementation kP(params);
    hardware::System system(hP, kP);
	physics::InterfacesHandlerImplementation interfacesHandler{params};
	physics::PrngParametersImplementation prngParameters{params};
	physics::PRNG prng{system, &prngParameters};
	
	//Operator for the test
	physics::fermionmatrix::MdagM_eo matrix(system, interfacesHandler.getInterface<physics::fermionmatrix::MdagM_eo>());
	//This configuration for the Ref.Code is the same as for example dks_input_5
	Gaugefield gf(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, std::string(SOURCEDIR) + "/ildg_io/conf.00200");
	
	hmc_float max,min;
	find_maxmin_eigenvalue(max, min, matrix, gf, system, interfacesHandler, 1.e-3, interfacesHandler.getAdditionalParameters<Staggeredfield_eo>());
	
	logger.info() << "mass squared = " << std::setprecision(16) << params.get_mass()*params.get_mass();
	logger.info() << " ref_max_eig = " << std::setprecision(16) << ref_max_eig;
	logger.info() << "     max_eig = " << std::setprecision(16) << max;
	logger.info() << " ref_min_eig = " << std::setprecision(16) << ref_min_eig;
	logger.info() << "     min_eig = " << std::setprecision(16) << min;
	
	BOOST_REQUIRE_SMALL(params.get_mass()*params.get_mass(), min);
	//The precision of this test is not so high because the method to calculate
	// the eigenvalue is not exactly the same as that implemented in the Ref.Code
	BOOST_CHECK_CLOSE(ref_max_eig, max, 1.e-3);
	BOOST_CHECK_CLOSE(ref_min_eig, min, 1.e-3);
}
