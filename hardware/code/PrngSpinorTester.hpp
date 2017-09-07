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
#include "SpinorStaggeredTester.hpp"
#include "prng.hpp"
#include "../../host_functionality/host_random.h" //@todo: remove this in the end!
#include "Kolmogorov_Smirnov.h"
#include "Normal_RNG_tests.h"

struct PrngSpinorTestParameters : public SpinorTestParameters
{
	PrngSpinorTestParameters(const LatticeExtents lE, const int iterationsIn) :
		TestParameters(lE, 10e-2), SpinorTestParameters(lE), iterations(iterationsIn) {}; // In calling the TestParameters ctor, the testPrecision is set to 10e-2, so as related tests can pass with a reasonable number of iterations!
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






struct PrngSpinorStaggeredTestParameters: public SpinorStaggeredTestParameters
{
	PrngSpinorStaggeredTestParameters(const LatticeExtents latticeExtentsIn, const int iterationsIn, const bool runKolmogorovSmirnov = true) :
		TestParameters(latticeExtentsIn, 10e-3), SpinorStaggeredTestParameters(latticeExtentsIn), iterations(iterationsIn),// In calling the TestParameters ctor, the testPrecision is set to 10e-3, so as related tests can pass with a reasonable number of iterations!
		runKolmogorovSmirnov(runKolmogorovSmirnov){};

	const unsigned int iterations;
	const bool runKolmogorovSmirnov;
};

struct PrngSpinorStaggeredTester: public SpinorStaggeredTester
{
	PrngSpinorStaggeredTester(const std::string kernelName, const ParameterCollection parameterCollection, const PrngSpinorStaggeredTestParameters & testParameters, const int numberOfElements, const ReferenceValues & referenceValues):
		SpinorStaggeredTester(kernelName, parameterCollection, testParameters, referenceValues),
		numberOfElements(numberOfElements), mean(0.), variance(0.),
		hostOutput(std::vector<su3vec>(numberOfElements * testParameters.iterations)),	testParameters(testParameters),
		hostSeed( parameterCollection.kernelParameters.getHostSeed() ),
		useSameRandomNumbers(parameterCollection.hardwareParameters.useSameRandomNumbers())
	{
		prng_init(hostSeed);
		prngStates = new hardware::buffers::PRNGBuffer(device, useSameRandomNumbers );
		auto codePrng = device->getPrngCode();
		codePrng->initialize(prngStates, hostSeed);
	}
	~PrngSpinorStaggeredTester()
	{
		calculateMean();
		calculateVariance();
		if (testParameters.runKolmogorovSmirnov)
			KolmogorovSmirnov();

		kernelResult.at(0) = mean;
		kernelResult.at(1) = sqrt(variance);

		delete prngStates;
	}
	double normalize(double valueIn)
	{
		return valueIn/= testParameters.iterations * numberOfElements * 6;
	}

	void calculateMean()
	{
		for (unsigned int i = 0; i < testParameters.iterations; i++) {
				if(i%100==0) logger.info() << "Run kernel for the " << i << "th time";
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

	void KolmogorovSmirnov()
	{
		vector<vector<hmc_float>> samples;
		vector<hmc_float> tmp;
		vector<hmc_float> tmp2;
		for(unsigned int i=0; i<testParameters.iterations; i++){
		  vector<hmc_float> tmp;
		  for(int j=0; j<numberOfElements; j++){
		    tmp2=reals_from_su3vec(hostOutput[i*numberOfElements+j]);
		    tmp.insert(tmp.end(),tmp2.begin(),tmp2.end());
		    tmp2.clear();
		  }
		  samples.push_back(tmp);
		  tmp.clear();
		}
		logger.info() << "Running Kolmogorov_Smirnov test (it should take half a minute)...";
		logger.info() << "Kolmogorov_Smirnov frequency (of K+): " << std::setprecision(16) << Kolmogorov_Smirnov(samples,0.,sqrt(0.5)) << " ---> It should be close 0.98";

		//Let us use the same sets of samples to make the mean test to 1,2 and 3 sigma
		//Note that in the test BOOST_CHECK is used.
		mean_test_multiple_set(samples,2.,0.,sqrt(0.5));
		mean_test_multiple_set(samples,3.,0.,sqrt(0.5));
		mean_test_multiple_set(samples,4.,0.,sqrt(0.5));
		//Let us use the same sets of samples to make the variance test to 1,2 and 3 sigma
		//Note that in the test BOOST_CHECK is used.
		variance_test_multiple_set(samples,2.,sqrt(0.5));
		variance_test_multiple_set(samples,3.,sqrt(0.5));
		variance_test_multiple_set(samples,4.,sqrt(0.5));
	}
protected:
	const int numberOfElements;
	double mean, variance;
	std::vector<su3vec> hostOutput;
	const PrngSpinorStaggeredTestParameters & testParameters;
	const hardware::buffers::PRNGBuffer* prngStates;
private:
	uint32_t hostSeed;
	bool useSameRandomNumbers;
};
