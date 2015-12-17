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

KernelTester::KernelTester (std::string kernelNameIn, const hardware::HardwareParametersInterface& hardwareParameters,
		const hardware::code::OpenClKernelParametersInterface& kernelParameters, const TestParameters testParams, const ReferenceValues rV) :
			testParameters(testParams),
			kernelResult(rV.size(),0),
			referenceValues(rV),
			hardwareParameters(&hardwareParameters),
			kernelParameters(&kernelParameters),
			kernelBuilder(nullptr)
{
	printKernelInformation(kernelNameIn);
	kernelBuilder = new hardware::OpenClCodeMockup( kernelParameters );
	system = new hardware::System(hardwareParameters, kernelParameters, *kernelBuilder);
	device = system->get_devices().at(0);
}

#include <boost/test/floating_point_comparison.hpp>
KernelTester::~KernelTester()
{
	//NOTE: Using "require" in boost throws an exception here, which should not happen in a destructor.
	for (int iteration = 0; iteration < (int) kernelResult.size(); iteration ++)
	{
		logger.info() << "compare result " << iteration;
		if (testParameters.typeOfComparison == ComparisonType::difference)
	    {
				logger.info() << std::setprecision(12) << "    Result = " << kernelResult[iteration];
				logger.info() << "Ref. Value = " << referenceValues[iteration];
				BOOST_CHECK_CLOSE(referenceValues[iteration], kernelResult[iteration], testParameters.testPrecision);
	    }
		else if (testParameters.typeOfComparison == ComparisonType::smallerThan)
	    {
				logger.info() << std::setprecision(12) << "    Result = " << kernelResult[iteration];
				logger.info() << "upper Bound = " << referenceValues[iteration];
				BOOST_CHECK_SMALL(kernelResult[iteration], referenceValues[iteration]);
	    }
	}

	if(system)
		delete system;
	if (kernelBuilder)
		delete kernelBuilder;
	device = nullptr;
}

