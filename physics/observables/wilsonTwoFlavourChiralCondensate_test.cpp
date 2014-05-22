/** @file
 * Unit test for the physics::observables::wilson::TwoFlavourChiralCondensate class
 *
 * Copyright 2014,Christopher Pinke
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

#include "wilsonTwoFlavourChiralCondensate.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics::observables::wilson::TwoFlavourChiralCondensate
#include <boost/test/unit_test.hpp>

#include "../../host_functionality/logger.hpp"
#include <stdexcept>

BOOST_AUTO_TEST_SUITE( BUILD )

	BOOST_AUTO_TEST_CASE( BUILD_1 )
	{
		const char * _params[] = {"foo", "--measure_pbp=true"};
		meta::Inputparameters params(2, _params);
		BOOST_REQUIRE_NO_THROW(physics::observables::wilson::TwoFlavourChiralCondensate tester(&params) );
	}
	
	void testLogicError(const char * _params[], int length )
	{
		meta::Inputparameters params(length, _params);
		BOOST_REQUIRE_THROW(physics::observables::wilson::TwoFlavourChiralCondensate tester(&params) , std::logic_error);
	}
	
	BOOST_AUTO_TEST_CASE( INV_ARGUMENT_1 )
	{
		const char * _params[] = {"foo", "--measure_pbp=false"};
		testLogicError(_params, 2);
	}
	
	BOOST_AUTO_TEST_CASE( INV_ARGUMENT_2 )
	{
		const char * _params[] = {"foo", "--measure_pbp=true", "--fermact=wilson", "--pbp_version=tm_one_end_trick"};
		testLogicError(_params, 4);
	}
	
	BOOST_AUTO_TEST_CASE( INV_ARGUMENT_3 )
	{
		const char * _params[] = {"foo", "--measure_pbp=true", "--fermact=clover", "--pbp_version=tm_one_end_trick"};
		testLogicError(_params, 4);
	}

	void testInvalidFermionAction(std::string actionName)
	{
		logger.info() << "Testing fermion action \"" + actionName + "\" for logic error";
		const char * standardParameters[] = {"foo", "--measure_pbp=true"};
		const char * commandLineParameters[] = {standardParameters[0], standardParameters[1], actionName.c_str() };
		
		meta::Inputparameters params(3, commandLineParameters);
		BOOST_REQUIRE_THROW(physics::observables::wilson::TwoFlavourChiralCondensate tester(&params) , std::logic_error);
	}
	
	BOOST_AUTO_TEST_CASE( INV_ARGUMENT_FERMION_ACTION )
	{
		std::vector<std::string> actionNames = {"clover", "tlsym", "iwasaki", "dbw2", "rooted_stagg"};
		for (int i = 0; i < (int) actionNames.size(); i++)
		{
			testInvalidFermionAction("--fermact=" + actionNames[i]);
		}
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( MEASURE )

	BOOST_AUTO_TEST_CASE( MEASURE_1 ) // equiv. to inverter test 29
	{
		const char * _params[] = {"foo", "--nt=4", "--ns=4", "--kappa=0.15", "--mu=4.", "--startcondition=cold", "--fermact=TWISTEDMASS", "--measure_pbp=true", "--num_sources=2", "--sourcetype=volume", "--sourcecontent=one", "--use_eo=false"};
		meta::Inputparameters params(11, _params);
		const hardware::System system(params);
		const physics::PRNG prng(system);
		const physics::lattices::Gaugefield gaugefield(system, prng);

		double result = physics::observables::wilson::measureChiralCondensate(&gaugefield, 0);
		double referenceValue = 4.86486486486488e-01;
		double testPrecision = 1e-8;
		BOOST_REQUIRE_CLOSE(result, referenceValue, testPrecision);
	}

BOOST_AUTO_TEST_SUITE_END()