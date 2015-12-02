/*
 * Copyright 2012, 2013, 2014 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
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

#include <iostream>
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE KernelTester_test
#include <boost/test/unit_test.hpp>

#include "mockups.hpp"

#include "kernelTester.hpp"

BOOST_AUTO_TEST_SUITE ( BUILD )

	BOOST_AUTO_TEST_CASE( BUILD_1 )
	{
		std::string nameOfKernel = "test";
		std::string nameOfInputfileThatExists = "kernelTesterEmpty_input";
		BOOST_CHECK_NO_THROW(KernelTester kernelTester(nameOfKernel, nameOfInputfileThatExists) );
	}

	BOOST_AUTO_TEST_CASE( BUILD_2 )
	{
		std::string nameOfKernel = "test";
		std::string nameOfInputfileThatDoesNotExist = "filethatdoesnotexist";
		BOOST_REQUIRE_THROW(KernelTester kernelTester(nameOfKernel, nameOfInputfileThatDoesNotExist), meta::Inputparameters::parse_aborted  );
	}

	BOOST_AUTO_TEST_CASE( INVALID_ARGUMENT_1 )
	{
		std::string nameOfKernel = "test";
		std::string nameOfInputfileThatExists = "kernelTesterEmpty_input";
		int maximumNumberOfReferenceValues = 2;
		BOOST_REQUIRE_THROW(KernelTester kernelTester(nameOfKernel, nameOfInputfileThatExists, maximumNumberOfReferenceValues + 1) , std::invalid_argument );
	}

	BOOST_AUTO_TEST_CASE( INVALID_ARGUMENT_2 )
	{
		std::string nameOfKernel = "test";
		std::string nameOfInputfileThatExists = "kernelTesterEmpty_input";
		int maximumNumberOfReferenceValues = 2;
		int minimalTypeOfComparision = 1;
		BOOST_REQUIRE_THROW(KernelTester kernelTester(nameOfKernel, nameOfInputfileThatExists, maximumNumberOfReferenceValues, minimalTypeOfComparision - 1), std::invalid_argument );
	}

	BOOST_AUTO_TEST_CASE( INVALID_ARGUMENT_3 )
	{
		std::string nameOfKernel = "test";
		std::string nameOfInputfileThatExists = "kernelTesterEmpty_input";
		int maximumNumberOfReferenceValues = 2;
		int maximumTypeOfComparision = 3;
		BOOST_REQUIRE_THROW(KernelTester kernelTester(nameOfKernel, nameOfInputfileThatExists, maximumNumberOfReferenceValues, maximumTypeOfComparision + 1), std::invalid_argument );
	}
	
	BOOST_AUTO_TEST_CASE( BUILD_3 )
	{
		const char * _params[] = {"foo"};
		meta::Inputparameters * parameter;
		hardware::System * system;
		hardware::Device * device;
		
		parameter = new meta::Inputparameters(1, _params);
		system = new hardware::System(*parameter);
		device = system->get_devices()[0];
		
		BOOST_CHECK_NO_THROW(KernelTester kernelTester(parameter, system, device) );
	}

	BOOST_AUTO_TEST_CASE( BUILD_4 )
	{
		const hardware::HardwareParametersMockup params(4,4);
		const hardware::code::OpenClKernelParametersMockup kernelParameters(4,4,true);
		const hardware::OpenClCodeMockup kernelBuilder(kernelParameters);
		const struct TestParameters testParams(LatticeExtents(4,4));
		BOOST_CHECK_NO_THROW(KernelTester kernelTester("testMockup", params, kernelParameters, testParams) );
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( DOUBLE )

	class TrivialKernelTester : public KernelTester {
	public:
		TrivialKernelTester(std::string kernelNameIn, std::string inputfileIn):
			KernelTester(kernelNameIn, inputfileIn) {
			kernelResult[0] = 1;
		}
	};

	BOOST_AUTO_TEST_CASE( TRIVIALKERNEL )
	{
		std::string nameOfKernel = "test";
		std::string nameOfInputfileThatExists = "kernelTester_input";
		TrivialKernelTester kernelTester(nameOfKernel, nameOfInputfileThatExists);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( DOUBLE_SMALL )

	class TrivialKernelTester : public KernelTester {
	public:
	  TrivialKernelTester(std::string kernelNameIn, std::string inputfileIn, double valueToCompare):
		  KernelTester(kernelNameIn, inputfileIn, 1, 2) {
			kernelResult[0] = valueToCompare;
		}
	};

BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES (TRIVIALKERNEL_1, 1)

	BOOST_AUTO_TEST_CASE( TRIVIALKERNEL_1 )
	{
		std::string nameOfKernel = "test";
		std::string nameOfInputfileThatExists = "kernelTester_input";
		double valueBiggerThanOne = 1.1;
		TrivialKernelTester kernelTester(nameOfKernel, nameOfInputfileThatExists, valueBiggerThanOne);
	}

	BOOST_AUTO_TEST_CASE( TRIVIALKERNEL_2 )
	{
		std::string nameOfKernel = "test";
		std::string nameOfInputfileThatExists = "kernelTester_input";
		double valueSmallerThanOne = 0.9;
		TrivialKernelTester kernelTester(nameOfKernel, nameOfInputfileThatExists, valueSmallerThanOne);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( COMPLEX )

	class TrivialKernelTester : public KernelTester {
	public:
		TrivialKernelTester(std::string kernelNameIn, std::string inputfileIn):
			KernelTester(kernelNameIn, inputfileIn, 2) {
			kernelResult[0] = 1.;
			kernelResult[1] = 2.;
		}
	};

	BOOST_AUTO_TEST_CASE( TRIVIALKERNEL )
	{
		std::string nameOfKernel = "test";
		std::string nameOfInputfileThatExists = "kernelTester_input";
		TrivialKernelTester kernelTester(nameOfKernel, nameOfInputfileThatExists);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( COMPLEX_ONE_REFERENCE_VALUE )

	class TrivialKernelTester : public KernelTester {
	public:
		TrivialKernelTester(std::string kernelNameIn, std::string inputfileIn):
			KernelTester(kernelNameIn, inputfileIn, 2, 3) {
			kernelResult[0] = 1.;
			kernelResult[1] = kernelResult[0];
		}
	};

	BOOST_AUTO_TEST_CASE( TRIVIALKERNEL )
	{
		std::string nameOfKernel = "test";
		std::string nameOfInputfileThatExists = "kernelTester_input";
		TrivialKernelTester kernelTester(nameOfKernel, nameOfInputfileThatExists);
	}

BOOST_AUTO_TEST_SUITE_END()
