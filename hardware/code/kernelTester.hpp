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

#include <boost/test/unit_test.hpp>
#include "../host_functionality/logger.hpp"
#include "testUtilities.hpp"
#include <vector>
#include "../system.hpp"
#include "../device.hpp"

class KernelTester {
public:
  KernelTester(std::string kernelNameIn, std::string inputfileIn, int numberOfValuesIn = 1, int typeOfComparision = 1);
	KernelTester(meta::Inputparameters * parameters, const hardware::System * system, hardware::Device * device);
	virtual ~KernelTester();	
	void setReferenceValuesToZero();
	
protected:
	double testPrecision;
	int typeOfComparision;
	std::vector<double> kernelResult;
	std::vector<double> referenceValue;
	bool allocatedObjects;
	
	meta::Inputparameters * parameters;
	const hardware::System * system;
	hardware::Device * device;
};

#endif /* KERNELTESTER_H_ */
