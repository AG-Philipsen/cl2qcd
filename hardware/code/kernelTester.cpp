/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
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

#include "kernelTester.hpp"
#include <boost/test/unit_test.hpp>

KernelTester::KernelTester(std::string kernelNameIn, std::string inputfileIn, int numberOfValuesIn, int typeOfComparisonIn):
  kernelResult(numberOfValuesIn, 0), referenceValue(numberOfValuesIn, 0), hardwareParameters(nullptr), kernelParameters(nullptr), kernelBuilder(nullptr)
{
	printKernelInformation(kernelNameIn);
	parameters = createParameters(inputfileIn).release();

	system = new hardware::System(*parameters);
	device = system->get_devices().at(0);
	allocatedObjects = true;

	testPrecision = parameters->get_solver_prec();

	for (int iteration = 0; iteration < (int) kernelResult.size(); iteration ++) {
		if(iteration == 0) {
			referenceValue[iteration] = parameters->get_test_ref_value();
		} else if(iteration == 1) {
			referenceValue[iteration] = parameters->get_test_ref_value2();
		} else {
			throw( std::invalid_argument("Can only set 2 reference values at the moment. Aborting...") );
		}
	}

	if ( (typeOfComparisonIn == 1) || (typeOfComparisonIn == 2)  || (typeOfComparisonIn == 3) )
	  {
	    typeOfComparison = typeOfComparisonIn;
	  } else
	  {
	    throw( std::invalid_argument("Do not recognize type of comparision. Aborting...") );
	  }
}

KernelTester::KernelTester(std::string kernelNameIn, std::vector<std::string> parameterStrings, size_t numberOfValuesIn, int typeOfComparisonIn, std::vector<double> expectedResult):
  kernelResult(numberOfValuesIn, 0), referenceValue(numberOfValuesIn, 0), hardwareParameters(nullptr), kernelParameters(nullptr), kernelBuilder(nullptr)
{
	printKernelInformation(kernelNameIn);
	parameters = createParameters(parameterStrings).release();

	system = new hardware::System(*parameters);
	device = system->get_devices()[0];
	allocatedObjects = true;
	
	testPrecision = parameters->get_solver_prec();

	if (expectedResult.size() == 0)
	{
		for (int iteration = 0; iteration < (int) kernelResult.size(); iteration ++) {
			if(iteration == 0) {
				referenceValue[iteration] = parameters->get_test_ref_value();
			} else if(iteration == 1) {
				referenceValue[iteration] = parameters->get_test_ref_value2();
			} else {
				throw( std::invalid_argument("Can only set 2 reference values at the moment. Aborting...") );
			}
		}
	}
	else
	{
		if( numberOfValuesIn != expectedResult.size() )
		{
			throw( std::invalid_argument("Number of arguments and size of expected results do not match. Aborting...") );
		}
		referenceValue = expectedResult;
	}

	if ( (typeOfComparisonIn == 1) || (typeOfComparisonIn == 2)  || (typeOfComparisonIn == 3) || (typeOfComparisonIn == 4) )
	  {
	    typeOfComparison = typeOfComparisonIn;
	  } else
	  {
	    throw( std::invalid_argument("Do not recognize type of comparision. Aborting...") );
	  }
}

KernelTester::KernelTester(meta::Inputparameters * parameters, const hardware::System * system, hardware::Device * device):
	testPrecision(1e-8), typeOfComparison(1), kernelResult(0, 0), referenceValue(0, 0), allocatedObjects(false), parameters(parameters), system(system), device(device), hardwareParameters(nullptr), kernelParameters(nullptr), kernelBuilder(nullptr)
{}

KernelTester::KernelTester (std::string kernelNameIn, const hardware::HardwareParametersInterface& hardwareParameters,
		const hardware::code::OpenClKernelParametersInterface& kernelParameters, struct TestParameters testParams) :
			kernelResult(testParams.numberOfValues,0),
			referenceValue(testParams.numberOfValues, 0),
			parameters(nullptr),
			hardwareParameters(&hardwareParameters),
			kernelParameters(&kernelParameters),
			kernelBuilder(nullptr)
{
	printKernelInformation(kernelNameIn);
	kernelBuilder = new hardware::OpenClCodeMockup( kernelParameters );
	system = new hardware::System(hardwareParameters, kernelParameters, *kernelBuilder);
	device = system->get_devices()[0];
	allocatedObjects = false;
	temporaryFlagForKernelTesterConstructorVersion = true;

	testPrecision = 10e-8; //todo: pass as arg

	if (testParams.referenceValue.size() == 0)
	{
		for (int iteration = 0; iteration < (int) kernelResult.size(); iteration ++) {
			if(iteration == 0) {
				referenceValue[iteration] = testParams.referenceValue.at(0);
			} else if(iteration == 1) {
				referenceValue[iteration] = testParams.referenceValue.at(1);
			} else {
				throw( std::invalid_argument("Can only set 2 reference values at the moment. Aborting...") );
			}
		}
	}
	else
	{
		if( testParams.numberOfValues != testParams.referenceValue.size() )
		{
			throw( std::invalid_argument("Number of arguments and size of expected results do not match. Aborting...") );
		}
		referenceValue = testParams.referenceValue;
	}

	//todo: the if and else can be removed if the enum is used anyway
	if ( (testParams.typeOfComparison == 1) || (testParams.typeOfComparison == 2)  || (testParams.typeOfComparison == 3) || (testParams.typeOfComparison == 4) )
	  {
	    typeOfComparison = testParams.typeOfComparison;
	  }
	else
	{
		logger.fatal() << testParams.typeOfComparison;
	    throw( std::invalid_argument("Do not recognise type of comparison. Aborting...") );
	}
}



#include <boost/test/floating_point_comparison.hpp>
KernelTester::~KernelTester()
{
  //NOTE: Using "require" in boost throws an exception here, which should not happen in a destructor.
	for (int iteration = 0; iteration < (int) kernelResult.size(); iteration ++) {
		logger.info() << "compare result " << iteration;
		if (typeOfComparison == ComparisonType::difference)
	    {
				logger.info() << std::setprecision(12) << "    Result = " << kernelResult[iteration];
				logger.info() << "Ref. Value = " << referenceValue[iteration];
				BOOST_CHECK_CLOSE(referenceValue[iteration], kernelResult[iteration], testPrecision);
	    }
		else if (typeOfComparison == ComparisonType::smallerThan)
	    {
				logger.info() << std::setprecision(12) << "    Result = " << kernelResult[iteration];
				logger.info() << "upper Bound = " << referenceValue[iteration];
				BOOST_CHECK_SMALL(kernelResult[iteration], referenceValue[iteration]);
	    }
		else if (typeOfComparison == ComparisonType::differenceToFirstReferenceValue)
	    {
				logger.info() << std::setprecision(12) << "    Result = " << kernelResult[iteration];
				logger.info() << "Ref. Value = " << referenceValue[0];
				BOOST_CHECK_CLOSE(referenceValue[0], kernelResult[iteration], testPrecision);
	    }
	}

	if(temporaryFlagForKernelTesterConstructorVersion)
	{
		delete system;
		delete kernelBuilder;
	}

	if(allocatedObjects)
	{
		delete parameters;
		delete system;
	}
	
	parameters = NULL;
	system = NULL;
	device = NULL;

}

void KernelTester::setReferenceValuesToZero()
{
	for (int iteration = 0; iteration < (int) referenceValue.size(); iteration ++) {
		referenceValue[iteration] = 0.;
	}
}
