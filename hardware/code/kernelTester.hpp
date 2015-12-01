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

#ifndef KERNELTESTER_H_
#define KERNELTESTER_H_

#include <vector>

#include <boost/test/unit_test.hpp>

#include "../../host_functionality/logger.hpp"
#include "testUtilities.hpp"
#include "../system.hpp"
#include "../device.hpp"
#include "mockups.hpp"

#include "latticeExtents.hpp"

enum ComparisonType{difference=1, smallerThan, differenceToFirstReferenceValue};

typedef std::vector<double> ReferenceValues;
static ReferenceValues defaultReferenceValues()
{
	return ReferenceValues{-1.23456};
}

struct TestParameters {
	std::vector<double> referenceValue;
	int ns;
	int nt;
	LatticeExtents latticeExtents;
	int typeOfComparison;
	size_t numberOfValues;

	//todo: introduce comparisonTypes here!
	TestParameters(std::vector<double> referenceValueIn, const LatticeExtents latticeExtentsIn):
		referenceValue(referenceValueIn), ns(latticeExtentsIn.ns), nt(latticeExtentsIn.nt), latticeExtents(latticeExtentsIn), typeOfComparison(ComparisonType::difference), numberOfValues(referenceValueIn.size()) {}
	TestParameters(std::vector<double> referenceValueIn, const LatticeExtents latticeExtentsIn, const int typeOfComparisonIn):
		referenceValue(referenceValueIn), ns(latticeExtentsIn.ns), nt(latticeExtentsIn.nt), latticeExtents(latticeExtentsIn),typeOfComparison(typeOfComparisonIn), numberOfValues(referenceValueIn.size()) {}
	TestParameters():
		referenceValue({0}), ns(4), nt(4), latticeExtents(LatticeExtents{4,4}),typeOfComparison(ComparisonType::difference), numberOfValues(referenceValue.size()) {}
};

struct ParameterCollection
{
	ParameterCollection(const hardware::HardwareParametersInterface & hardwareParametersIn, const hardware::code::OpenClKernelParametersInterface & kernelParametersIn):
		hardwareParameters(hardwareParametersIn), kernelParameters(kernelParametersIn) {};
	const hardware::HardwareParametersInterface & hardwareParameters;
	const hardware::code::OpenClKernelParametersInterface & kernelParameters;
};

class KernelTester {
public:
	KernelTester(std::string kernelNameIn, const hardware::HardwareParametersInterface&, const hardware::code::OpenClKernelParametersInterface&, struct TestParameters);
	KernelTester(std::string kernelNameIn, const hardware::HardwareParametersInterface&, const hardware::code::OpenClKernelParametersInterface&, struct TestParameters, const ReferenceValues);
	//todo: remove these constructors
	KernelTester(std::string kernelNameIn, std::string inputfileIn, int numberOfValuesIn = 1, int typeOfComparison = 1);
	KernelTester(std::string kernelNameIn, std::vector<std::string> parameterStrings, size_t numberOfValuesIn = 1, int typeOfComparison = 1, std::vector<double> result = std::vector<double>() );
	KernelTester(meta::Inputparameters * parameters, const hardware::System * system, hardware::Device * device);
	virtual ~KernelTester();
	void setReferenceValuesToZero();
	
protected:
	double testPrecision;
	int typeOfComparison; //todo: introduce comparisonTypes here!
	std::vector<double> kernelResult;
	std::vector<double> referenceValue; //@todo: must be ReferenceValues
	bool allocatedObjects; //todo: remove
	
	meta::Inputparameters * parameters; //todo: remove
	const hardware::System * system;
	hardware::Device * device;
	const hardware::HardwareParametersInterface * hardwareParameters;
	const hardware::code::OpenClKernelParametersInterface * kernelParameters;
	const hardware::OpenClCode * kernelBuilder; //todo: remove

private:
	bool temporaryFlagForKernelTesterConstructorVersion = false; //todo: remove
};

#endif /* KERNELTESTER_H_ */
