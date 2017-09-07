/** @file
 * Unit tests for the physics::observables::wilson::TwoFlavourChiralCondensate class
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
#include <boost/lexical_cast.hpp>

#include "../../interfaceImplementations/interfacesHandler.hpp"
#include "../../interfaceImplementations/hardwareParameters.hpp"
#include "../../interfaceImplementations/openClKernelParameters.hpp"

BOOST_AUTO_TEST_SUITE( BUILD )

	void testLogicError(const char * _params[], int length )
	{
		meta::Inputparameters params(length, _params);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::PrngParametersImplementation prngParameters{params};
		const physics::PRNG prng{system, &prngParameters};
		const physics::lattices::Gaugefield gaugefield(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng);

		
		BOOST_REQUIRE_THROW(physics::observables::wilson::measureTwoFlavourChiralCondensateAndWriteToFile(&gaugefield, 0, interfacesHandler) , std::logic_error);
	}
	
	BOOST_AUTO_TEST_CASE( INV_ARGUMENT_1 )
	{
		const char * _params[] = {"foo", "--measure_pbp=false"};
		testLogicError(_params, 2);
	}

	void testInvalidFermionActionAndVersion(std::string actionName, std::string version = "--pbp_version=std")
	{
		logger.info() << "Testing fermion action \"" + actionName + "\" and chiral condensate version \"" + version + "\" for logic error";
		const char * standardParameters[] = {"foo", "--measure_pbp=true"};
		const char * commandLineParameters[] = {standardParameters[0], standardParameters[1], actionName.c_str() , version.c_str()};
		
		meta::Inputparameters params(4, commandLineParameters);
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
		physics::PrngParametersImplementation prngParameters{params};
		const physics::PRNG prng{system, &prngParameters};
		const physics::lattices::Gaugefield gaugefield(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng);
		
		BOOST_REQUIRE_THROW(physics::observables::wilson::measureTwoFlavourChiralCondensateAndWriteToFile(&gaugefield, 0, interfacesHandler) , std::logic_error);
	}
	
	std::vector<std::string> actionNames = {"clover", "tlsym", "iwasaki", "dbw2", "rooted_stagg"};
	
	BOOST_AUTO_TEST_CASE( INV_ARGUMENT_FERMION_ACTION )
	{
		for (int i = 0; i < (int) actionNames.size(); i++)
		{
			testInvalidFermionActionAndVersion("--fermact=" + actionNames[i]);
		}
	}
	
	BOOST_AUTO_TEST_CASE( INV_ARGUMENT_CHIRAL_CONDENSATE_VERSION )
	{
		actionNames.push_back("wilson");
		for (int i = 0; i < (int) actionNames.size(); i++)
		{
			testInvalidFermionActionAndVersion("--fermact=" + actionNames[i], "--pbp_version=tm_one_end_trick");
		}
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( MEASURE )

void testMeasurement(std::vector<double> referenceValues, int numberOfSources, std::string fermactOptionIn, std::string sourceTypeOptionIn, std::string pbpVersionOptionIn, std::string eoOptionIn = "false", std::string startconditionOptionIn = "cold", std::string sourcefileOptionIn = "conf.00200")
{
		std::string numberOfSources_option = "--num_sources=" + boost::lexical_cast<std::string>( numberOfSources );
		std::string fermactOption = "--fermact=" + fermactOptionIn;
		std::string sourceTypeOption = "--sourcetype=" + sourceTypeOptionIn;
		std::string sourceContentOption = "--sourcecontent=one";
		std::string pbpVersionOption = "--pbp_version=" + pbpVersionOptionIn;
		std::string startconditionOption = "--startcondition=" + startconditionOptionIn;
		std::string sourcefileOption = "--sourcefile=" + sourcefileOptionIn;
		std::string eoOption = "--use_eo=" + eoOptionIn;

		int numberOfOptions= 14;

		int numberOfParametersFromBoost =       boost::unit_test::framework::master_test_suite().argc;
		std::vector<std::string> boostOptions;
		if (numberOfParametersFromBoost > 1)
		  {
		    for (int i = 1; i < numberOfParametersFromBoost; i++)
		      {
			boostOptions.push_back( boost::unit_test::framework::master_test_suite().argv[i]);
		      }
		    numberOfOptions ++;
		  }
		else
		  {
		    boostOptions.push_back("");
		  }

		const char * _params[] = {"foo", "--nt=4", "--ns=4", "--kappa=0.15", "--mu=4.", "--measure_pbp=true", fermactOption.c_str(), sourceTypeOption.c_str(), sourceContentOption.c_str(),  numberOfSources_option.c_str(), pbpVersionOption.c_str(), eoOption.c_str(), startconditionOption.c_str(), sourcefileOption.c_str(), boostOptions[0].c_str()};

		meta::Inputparameters params(numberOfOptions, _params);
	    hardware::HardwareParametersImplementation hP(&params);
	    hardware::code::OpenClKernelParametersImplementation kP(params);
	    hardware::System system(hP, kP);
		physics::InterfacesHandlerImplementation interfacesHandler{params};
		physics::PrngParametersImplementation prngParameters{params};
		const physics::PRNG prng{system, &prngParameters};
		const physics::lattices::Gaugefield gaugefield(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng);

		std::vector<double> results;
		results = physics::observables::wilson::measureTwoFlavourChiralCondensateAndWriteToFile(&gaugefield, "conf.test", interfacesHandler);
		
		BOOST_REQUIRE_EQUAL(numberOfSources, (int) results.size() );

		double testPrecision = 1e-8;		
		for (int i = 0; i < (int) referenceValues.size(); i++)
		{
			BOOST_REQUIRE_CLOSE(results[i], referenceValues[i], testPrecision);
		}
}

        BOOST_AUTO_TEST_CASE( MEASURE_NONEO )
	{
		int numberOfSources = 9;
		std::vector<double> referenceValues(numberOfSources, 4.86486486486488e-01);
	
		testMeasurement(referenceValues, numberOfSources, (std::string) "twistedmass", (std::string) "volume", (std::string) "std", (std::string) "false");
	}

        BOOST_AUTO_TEST_CASE( MEASURE_EO )
	{
		int numberOfSources = 9;
		std::vector<double> referenceValues(numberOfSources, 4.86486486486488e-01);
	
		testMeasurement(referenceValues, numberOfSources, (std::string) "twistedmass", (std::string) "volume", (std::string) "std", (std::string) "true");
	}

        BOOST_AUTO_TEST_CASE( MEASURE_ONE_END_TRICK )
	{
		int numberOfSources = 9;
		std::vector<double> referenceValues(numberOfSources, 4.86486486486488e-01);
	
		testMeasurement(referenceValues, numberOfSources, (std::string) "twistedmass", (std::string) "volume", (std::string) "tm_one_end_trick");
	}

BOOST_AUTO_TEST_CASE( MEASURE_Z_SLICE )
	{
		int numberOfSources = 9;
		std::vector<double> referenceValues(numberOfSources, 1.16971963846964e-01);
	
		testMeasurement(referenceValues, numberOfSources, (std::string) "twistedmass", (std::string) "zslice", (std::string) "std");
	}

BOOST_AUTO_TEST_CASE( MEASURE_CONFIGURATION_1 )
	{
		int numberOfSources = 9;
		std::vector<double> referenceValues(numberOfSources, 6.86307941352032e-02);
	
		testMeasurement(referenceValues, numberOfSources, (std::string) "twistedmass", (std::string) "zslice", (std::string) "std", (std::string) "false", (std::string) "continue");
	}

BOOST_AUTO_TEST_CASE( MEASURE_CONFIGURATION_2 )
	{
		int numberOfSources = 9;
		std::vector<double> referenceValues(numberOfSources, 1.35706953188924e-01);
	
		testMeasurement(referenceValues, numberOfSources, (std::string) "wilson", (std::string) "zslice", (std::string) "std", (std::string) "false", (std::string) "continue");
	}

BOOST_AUTO_TEST_SUITE_END()
