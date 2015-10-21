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

typedef std::vector<double> referenceValues;
struct TestParameters {
	std::vector<double> referenceValue;
	int ns;
	int nt;
	int typeOfComparison;
	size_t numberOfValues;

	TestParameters(std::vector<double> referenceValueIn, int nsIn, int ntIn): referenceValue(referenceValueIn), ns(nsIn), nt(ntIn), typeOfComparison(1), numberOfValues(referenceValueIn.size()) {}
	TestParameters(): referenceValue({-1234}), ns(4), nt(4), typeOfComparison(1), numberOfValues(referenceValue.size()) {}
};

class KernelTester {
public:
	KernelTester(std::string kernelNameIn, std::string inputfileIn, int numberOfValuesIn = 1, int typeOfComparison = 1);
	KernelTester(std::string kernelNameIn, std::vector<std::string> parameterStrings, size_t numberOfValuesIn = 1, int typeOfComparison = 1, std::vector<double> result = std::vector<double>() );
	KernelTester(meta::Inputparameters * parameters, const hardware::System * system, hardware::Device * device);
	KernelTester(std::string kernelNameIn, const hardware::HardwareParametersInterface&, const hardware::code::OpenClKernelParametersInterface&,
			struct TestParameters);
	virtual ~KernelTester();
	void setReferenceValuesToZero();
	
protected:
	double testPrecision;
	int typeOfComparison;
	std::vector<double> kernelResult;
	std::vector<double> referenceValue;
	bool allocatedObjects;
	
	meta::Inputparameters * parameters;
	const hardware::System * system;
	hardware::Device * device;
	const hardware::HardwareParametersInterface * hardwareParameters;
	const hardware::code::OpenClKernelParametersInterface * kernelParameters;
	const hardware::OpenClCode * kernelBuilder;

private:
	bool temporaryFlagForKernelTesterConstructorVersion = false;
};

#endif /* KERNELTESTER_H_ */
