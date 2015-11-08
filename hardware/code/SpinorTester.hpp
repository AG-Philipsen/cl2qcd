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
#include "../../physics/prng.hpp" //todo: remove
#include "spinors.hpp"
#include "complex.hpp"

enum SpinorFillType{ zero, one, zeroOne, oneZero, ascendingReal, ascendingComplex};
typedef std::vector<SpinorFillType> SpinorFillTypes;
typedef std::vector<hmc_complex> ComplexNumbers;
typedef size_t NumberOfSpinors;

struct LatticeExtends
{
	LatticeExtends() : ns(4), nt(4) {};
	LatticeExtends(const int nsIn, const int ntIn): ns(nsIn), nt(ntIn) {};
	LatticeExtends(const unsigned int nsIn, const unsigned int ntIn): ns(nsIn), nt(ntIn) {};
	const unsigned int ns;
	const unsigned int nt;
};

static int calculateLatticeVolume(const int nsIn, const int ntIn) noexcept
{
	return 	nsIn * nsIn * nsIn * ntIn;
}

static int calculateSpinorfieldSize(const int nsIn, const int ntIn) noexcept
{
	return 	calculateLatticeVolume(nsIn, ntIn);
}

static int calculateSpinorfieldSize(const LatticeExtends latticeExtendsIn) noexcept
{
	return 	calculateLatticeVolume(latticeExtendsIn.ns, latticeExtendsIn.nt);
}

static int calculateEvenOddSpinorfieldSize(const int nsIn, const int ntIn) noexcept
{
	return 	calculateSpinorfieldSize(nsIn, ntIn) / 2;
}

static int calculateEvenOddSpinorfieldSize(const LatticeExtends latticeExtendsIn) noexcept
{
	return 	calculateSpinorfieldSize(latticeExtendsIn.ns, latticeExtendsIn.nt) / 2;
}

struct SpinorTestParameters: public TestParameters
{
	SpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtends latticeExtendsIn, const SpinorFillTypes fillTypesIn, const bool needEvenOddIn) :
		TestParameters(referenceValuesIn, latticeExtendsIn.ns, latticeExtendsIn.nt), needEvenOdd(needEvenOddIn), fillTypes(fillTypesIn) {};
	SpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtends latticeExtendsIn) :
		TestParameters(referenceValuesIn, latticeExtendsIn.ns, latticeExtendsIn.nt), needEvenOdd(false), fillTypes(SpinorFillType::one) {};
	SpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtends latticeExtendsIn, const bool needEvenOddIn) :
		TestParameters(referenceValuesIn, latticeExtendsIn.ns, latticeExtendsIn.nt), needEvenOdd(needEvenOddIn), fillTypes(SpinorFillType::one) {};
	SpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtends latticeExtendsIn, const SpinorFillTypes fillTypesIn, const bool needEvenOddIn, const int typeOfComparisionIn) :
		TestParameters(referenceValuesIn, latticeExtendsIn.ns, latticeExtendsIn.nt, typeOfComparisionIn), needEvenOdd(needEvenOddIn), fillTypes(fillTypesIn) {};
	SpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtends latticeExtendsIn, const int typeOfComparisionIn, const bool needEvenOddIn) :
		TestParameters(referenceValuesIn, latticeExtendsIn.ns, latticeExtendsIn.nt, typeOfComparisionIn), needEvenOdd(needEvenOddIn), fillTypes(SpinorFillType::one) {};
	SpinorTestParameters() : TestParameters(), needEvenOdd(false) {};

	int getSpinorfieldSize() const { return calculateSpinorfieldSize(ns, nt); } ;
	int getEvenOddSpinorfieldSize() const { return calculateEvenOddSpinorfieldSize(ns, nt); } ;
	int getSpinorfieldSize(const LatticeExtends latticeExtendsIn) const { return calculateSpinorfieldSize(latticeExtendsIn); } ;
	int getEvenOddSpinorfieldSize(const LatticeExtends latticeExtendsIn) const { return calculateEvenOddSpinorfieldSize(latticeExtendsIn); } ;

	const bool needEvenOdd;
	const SpinorFillTypes fillTypes;
};

struct NonEvenOddSpinorTestParameters : public SpinorTestParameters
{
	NonEvenOddSpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtends latticeExtendsIn, const SpinorFillTypes fillTypesIn) :
		SpinorTestParameters(referenceValuesIn, latticeExtendsIn, fillTypesIn, false) {};
	NonEvenOddSpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtends latticeExtendsIn) :
		SpinorTestParameters(referenceValuesIn,latticeExtendsIn, false) {};
	NonEvenOddSpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtends latticeExtendsIn, const SpinorFillTypes fillTypesIn, const ComparisonType typeOfComparisionIn) :
		SpinorTestParameters(referenceValuesIn, latticeExtendsIn, fillTypesIn, false, typeOfComparisionIn) {};
	NonEvenOddSpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtends latticeExtendsIn, const ComparisonType typeOfComparisionIn) :
		SpinorTestParameters(referenceValuesIn, latticeExtendsIn, SpinorFillTypes{SpinorFillType::one}, false, typeOfComparisionIn) {};
	NonEvenOddSpinorTestParameters() : SpinorTestParameters() {};
};

struct EvenOddSpinorTestParameters : public SpinorTestParameters
{
	EvenOddSpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtends latticeExtendsIn, const SpinorFillTypes fillTypesIn) :
		SpinorTestParameters(referenceValuesIn, latticeExtendsIn, fillTypesIn, true) {};
	EvenOddSpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtends latticeExtendsIn) :
		SpinorTestParameters(referenceValuesIn,latticeExtendsIn, true) {};
	EvenOddSpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtends latticeExtendsIn, const SpinorFillTypes fillTypesIn, const ComparisonType typeOfComparisionIn) :
		SpinorTestParameters(referenceValuesIn, latticeExtendsIn, fillTypesIn, true, typeOfComparisionIn) {};
	EvenOddSpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtends latticeExtendsIn, const ComparisonType typeOfComparisionIn) :
		SpinorTestParameters(referenceValuesIn, latticeExtendsIn, SpinorFillTypes{SpinorFillType::one}, true, typeOfComparisionIn) {};
	EvenOddSpinorTestParameters() : SpinorTestParameters() {};
};

class SpinorTester : public KernelTester {
public:
	SpinorTester(std::string kernelName, const ParameterCollection,	const SpinorTestParameters & testParameters );
	//todo: remove these constructors
	SpinorTester(std::string kernelName, std::string inputfileIn, int numberOfValues = 1, int typeOfComparision = 1);
	SpinorTester(std::string kernelName,  std::vector<std::string> parameterStrings, int numberOfValues = 1, int typeOfComparision = 1, std::vector<double> expectedResult = std::vector<double> ());
	SpinorTester(meta::Inputparameters * parameters, const hardware::System * system, hardware::Device * device);
	~SpinorTester();
	
protected:
	std::string getSpecificInputfile(std::string inputfileIn); //todo: remove
	
	bool allocatedObjects; //todo: remove

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
	spinor * createSpinorfieldWithOnesAndZerosDependingOnSiteParity(const bool fillEvenSites);
	spinor * createSpinorfieldWithOnesAndMinusOneForGamma5Use(size_t numberOfElements);	
	void fillTwoSpinorBuffersDependingOnParity(const hardware::buffers::Spinor * in1, const hardware::buffers::Spinor * in2);
	void fillTwoSpinorfieldsDependingOnParity(spinor * sf_in1, spinor * sf_in2, int size);
	void fill_with_one_eo(spinor * in, const int, const bool);
	hmc_float count_sf(spinor * in, int size);
	hmc_float calc_var(hmc_float in, hmc_float mean);
	hmc_float calc_var_sf(spinor * in, int size, hmc_float sum);
	void calcSquarenormAndStoreAsKernelResult(const hardware::buffers::Plain<spinor> * in);
	void calcSquarenormEvenOddAndStoreAsKernelResult(const hardware::buffers::Spinor * in);
	void fillTwoSpinorfieldsWithRandomNumbers(spinor * sf_in1, spinor * sf_in2, int size, int seed = 123456);
	
	void setMembers(); //todo: remove
	void setMembersNew();
	
	const hardware::code::Spinors * code;
	//todo: remove
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
