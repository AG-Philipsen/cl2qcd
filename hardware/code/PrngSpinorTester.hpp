/*
 * Copyright 2015 Christopher Pinke
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

#pragma once

#include "SpinorTester.hpp"
#include "prng.hpp"
#include "../../host_functionality/host_random.h" //@todo: remove this in the end!

struct PrngSpinorTestParameters : public SpinorTestParameters
{
	PrngSpinorTestParameters(const LatticeExtents lE) :
		TestParameters(lE, ComparisonType::smallerThan), SpinorTestParameters(lE, ComparisonType::smallerThan), iterations(100) {};
	PrngSpinorTestParameters(const LatticeExtents lE, const int iterationsIn) :
		TestParameters(lE, ComparisonType::smallerThan), SpinorTestParameters(lE, ComparisonType::smallerThan), iterations(iterationsIn) {};
	const unsigned int iterations;
};

struct PrngSpinorTester: public SpinorTester
{
	PrngSpinorTester(const std::string kernelName, const ParameterCollection parameterCollection, const PrngSpinorTestParameters & testParameters, const int numberOfElements, const ReferenceValues rV):
				SpinorTester(kernelName, parameterCollection, testParameters, rV),
				numberOfElements(numberOfElements), mean(0.), variance(0.),
				hostOutput(std::vector<spinor> (numberOfElements * testParameters.iterations)),	testParameters(testParameters),
				hostSeed( parameterCollection.kernelParameters.getHostSeed() ),
				useSameRandomNumbers(parameterCollection.hardwareParameters.useSameRandomNumbers())
	{
		prng_init(hostSeed);
		prngStates = new hardware::buffers::PRNGBuffer(device, useSameRandomNumbers );
		auto codePrng = device->getPrngCode();
		codePrng->initialize(prngStates, hostSeed);
	}

	~PrngSpinorTester()
	{
		calculateMean();
		calculateVariance();

		kernelResult.at(0) = mean;
		kernelResult.at(1) = sqrt(variance);

		delete prngStates;
	}

	double normalize(double valueIn)
	{
		return valueIn/= testParameters.iterations * numberOfElements * 24;
	}

	void calculateMean()
	{
		for (unsigned int i = 0; i < testParameters.iterations; i++) {
			mean += count_sf(&hostOutput[i * numberOfElements], numberOfElements);
		}
		mean = normalize(mean);
	}

	void calculateVariance()
	{
		for (unsigned int i = 0; i < testParameters.iterations; i++) {
			variance += calc_var_sf(&hostOutput[i * numberOfElements], numberOfElements, mean);
		}
		variance = normalize(variance);
	}
protected:
	const int numberOfElements;
	double mean, variance;
	std::vector<spinor> hostOutput;
	const PrngSpinorTestParameters & testParameters;
	const hardware::buffers::PRNGBuffer* prngStates;
private:
	uint32_t hostSeed;
	bool useSameRandomNumbers;
};

