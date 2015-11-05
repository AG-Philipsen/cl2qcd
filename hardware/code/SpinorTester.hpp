/*
 * Copyright 2014 Christopher Pinke
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

#ifndef SPINORTESTER_HPP_
#define SPINORTESTER_HPP_

#include "kernelTester.hpp"

#include "../../meta/util.hpp"
#include "../../host_functionality/host_random.h"
#include "../../physics/prng.hpp"
#include "spinors.hpp"
#include "complex.hpp"

enum SpinorFillType{ zero, one, zeroOne, oneZero, ascendingReal, ascendingComplex};
typedef std::vector<SpinorFillType> SpinorFillTypes;
typedef std::vector<hmc_complex> ComplexNumbers;
typedef int numberOfSpinors;

struct SpinorTestParameters: public TestParameters
{
	SpinorTestParameters(const ReferenceValues referenceValuesIn, const int nsIn, const int ntIn, const SpinorFillTypes fillTypesIn, const bool isEvenOddIn) :
		TestParameters(referenceValuesIn, nsIn, ntIn), isEvenOdd(isEvenOddIn), fillTypes(fillTypesIn), typeOfComparison(1) {};
	SpinorTestParameters(const ReferenceValues referenceValuesIn, const int nsIn, const int ntIn, const int typeOfComparisonIn) :
		TestParameters(referenceValuesIn, nsIn, ntIn), isEvenOdd(false), fillTypes(SpinorFillType::one), typeOfComparison(typeOfComparisonIn) {};
	SpinorTestParameters(const ReferenceValues referenceValuesIn, const int nsIn, const int ntIn, const bool isEvenOddIn) :
		TestParameters(referenceValuesIn, nsIn, ntIn), isEvenOdd(isEvenOddIn), fillTypes(SpinorFillType::one), typeOfComparison(1) {};
	const bool isEvenOdd;
	const SpinorFillTypes fillTypes;
	const int typeOfComparison;
};

class SpinorTester : public KernelTester {
public:
	SpinorTester(std::string kernelName, std::string inputfileIn, int numberOfValues = 1, int typeOfComparision = 1);
	SpinorTester(std::string kernelName,  std::vector<std::string> parameterStrings, int numberOfValues = 1, int typeOfComparision = 1, std::vector<double> expectedResult = std::vector<double> ());
	SpinorTester(meta::Inputparameters * parameters, const hardware::System * system, hardware::Device * device);
	SpinorTester(std::string kernelName, const hardware::HardwareParametersInterface &, const hardware::code::OpenClKernelParametersInterface &,
			const SpinorTestParameters & testParameters );
	~SpinorTester();
	
protected:
	std::string getSpecificInputfile(std::string inputfileIn);
	
	bool allocatedObjects;

	spinor * createSpinorfield( SpinorFillType );
	spinor * createSpinorfield(size_t numberOfElements, int seed = 123456);
	void fillTwoSpinorBuffers(const hardware::buffers::Spinor * in1, const hardware::buffers::Spinor * in2, int seed = 123456);
	void fill_with_one(spinor * in, int size);
	void fill_with_zero_one(spinor * in, int size);
	void fill_with_one_zero(spinor * in, int size);
	void fill_with_ascending(spinor * in, int size);
	void fillWithAscendingComplex(spinor * in, int size);
	void fill_with_one_minusone_for_gamma5_use(spinor * in, int size);
	void fill_with_random(spinor * in, int size, int seed);
	spinor * createSpinorfieldWithOnesAndZerosDependingOnSiteParity();
	spinor * createSpinorfieldWithOnesAndMinusOneForGamma5Use(size_t numberOfElements);	
	void fillTwoSpinorBuffersDependingOnParity(const hardware::buffers::Spinor * in1, const hardware::buffers::Spinor * in2);
	void fillTwoSpinorfieldsDependingOnParity(spinor * sf_in1, spinor * sf_in2, int size);
	void fill_with_one_eo(spinor * in, int size, bool eo);
	hmc_float count_sf(spinor * in, int size);
	hmc_float calc_var(hmc_float in, hmc_float mean);
	hmc_float calc_var_sf(spinor * in, int size, hmc_float sum);
	void calcSquarenormAndStoreAsKernelResult(const hardware::buffers::Plain<spinor> * in);
	void calcSquarenormEvenOddAndStoreAsKernelResult(const hardware::buffers::Spinor * in);
	void fillTwoSpinorfieldsWithRandomNumbers(spinor * sf_in1, spinor * sf_in2, int size, int seed = 123456);
	
	void setMembers();
	void setMembersNew();
	
	const hardware::code::Spinors * code;
	const physics::ParametersPrng_fromMetaInputparameters prngParameters;
	physics::PRNG * prng;

	hardware::buffers::Plain<double> * doubleBuffer;
	
	//todo: most of these must go away!
	size_t spinorfieldElements;
	size_t spinorfieldEvenOddElements;
	bool useRandom;
	bool evenOrOdd;
	bool calcVariance;
	hmc_complex alpha_host;
	hmc_complex beta_host;
	int iterations;
};

#endif
