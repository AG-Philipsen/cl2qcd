/*
 * Copyright 2014, 2015 Christopher Pinke
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

#include "kernelTester.hpp"
#include "spinors.hpp"

//@todo: check if in the end, ascending is used at all. If not, remove and rename ascendingComplex to ascending
enum SpinorFillType{ zero, one, zeroOne, oneZero, ascendingReal, ascendingComplex};
typedef std::vector<SpinorFillType> SpinorFillTypes;
typedef std::vector<hmc_complex> ComplexNumbers;
typedef size_t NumberOfSpinors;

static int calculateSpinorfieldSize(const int nsIn, const int ntIn) noexcept
{
	return 	calculateLatticeVolume(nsIn, ntIn);
}

static int calculateSpinorfieldSize(const LatticeExtents latticeExtendsIn) noexcept
{
	return 	calculateLatticeVolume(latticeExtendsIn.ns, latticeExtendsIn.nt);
}

static int calculateEvenOddSpinorfieldSize(const int nsIn, const int ntIn) noexcept
{
	return 	calculateSpinorfieldSize(nsIn, ntIn) / 2;
}

static int calculateEvenOddSpinorfieldSize(const LatticeExtents latticeExtendsIn) noexcept
{
	return 	calculateSpinorfieldSize(latticeExtendsIn.ns, latticeExtendsIn.nt) / 2;
}

struct SpinorTestParameters: public TestParameters
{
	SpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn, const bool needEvenOddIn) :
		TestParameters(referenceValuesIn, latticeExtendsIn), needEvenOdd(needEvenOddIn), fillTypes(fillTypesIn) {};
	SpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn) :
		TestParameters(referenceValuesIn, latticeExtendsIn), needEvenOdd(false), fillTypes(SpinorFillType::one) {};
	SpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn, const bool needEvenOddIn) :
		TestParameters(referenceValuesIn, latticeExtendsIn), needEvenOdd(needEvenOddIn), fillTypes(SpinorFillType::one) {};
	SpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn, const bool needEvenOddIn, const int typeOfComparisionIn) :
		TestParameters(referenceValuesIn, latticeExtendsIn, typeOfComparisionIn), needEvenOdd(needEvenOddIn), fillTypes(fillTypesIn) {};
	SpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn, const int typeOfComparisionIn, const bool needEvenOddIn) :
		TestParameters(referenceValuesIn, latticeExtendsIn, typeOfComparisionIn), needEvenOdd(needEvenOddIn), fillTypes(SpinorFillType::one) {};
	SpinorTestParameters() : TestParameters(), needEvenOdd(false) {};

	int getSpinorfieldSize() const { return calculateSpinorfieldSize(ns, nt); } ; //@todo: remove?
	int getEvenOddSpinorfieldSize() const { return calculateEvenOddSpinorfieldSize(ns, nt); } ; //@todo: remove?
	int getSpinorfieldSize(const LatticeExtents latticeExtendsIn) const { return calculateSpinorfieldSize(latticeExtendsIn); } ; //@todo: remove
	int getEvenOddSpinorfieldSize(const LatticeExtents latticeExtendsIn) const { return calculateEvenOddSpinorfieldSize(latticeExtendsIn); } ; //@todo: remove

	const bool needEvenOdd; //@todo: remove, this should be in the fields themselves!
	const SpinorFillTypes fillTypes;
};

struct NonEvenOddSpinorTestParameters : public SpinorTestParameters
{
	NonEvenOddSpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn) :
		SpinorTestParameters(referenceValuesIn, latticeExtendsIn, fillTypesIn, false) {};
	NonEvenOddSpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn, const SpinorFillType fillTypeIn) :
		SpinorTestParameters(referenceValuesIn, latticeExtendsIn, SpinorFillTypes{fillTypeIn}, false) {};
	NonEvenOddSpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn) :
		SpinorTestParameters(referenceValuesIn,latticeExtendsIn, false) {};
	NonEvenOddSpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn, const ComparisonType typeOfComparisionIn) :
		SpinorTestParameters(referenceValuesIn, latticeExtendsIn, fillTypesIn, false, typeOfComparisionIn) {};
	NonEvenOddSpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn, const ComparisonType typeOfComparisionIn) :
		SpinorTestParameters(referenceValuesIn, latticeExtendsIn, SpinorFillTypes{SpinorFillType::one}, false, typeOfComparisionIn) {};
	NonEvenOddSpinorTestParameters() : SpinorTestParameters() {};
};

struct EvenOddSpinorTestParameters : public SpinorTestParameters
{
	EvenOddSpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn) :
		SpinorTestParameters(referenceValuesIn, latticeExtendsIn, fillTypesIn, true) {};
	EvenOddSpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn, const SpinorFillType fillTypeIn) :
		SpinorTestParameters(referenceValuesIn, latticeExtendsIn, SpinorFillTypes{fillTypeIn}, true) {};
	EvenOddSpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn) :
		SpinorTestParameters(referenceValuesIn,latticeExtendsIn, true) {};
	EvenOddSpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn, const ComparisonType typeOfComparisionIn) :
		SpinorTestParameters(referenceValuesIn, latticeExtendsIn, fillTypesIn, true, typeOfComparisionIn) {};
	EvenOddSpinorTestParameters(const ReferenceValues referenceValuesIn, const LatticeExtents latticeExtendsIn, const ComparisonType typeOfComparisionIn) :
		SpinorTestParameters(referenceValuesIn, latticeExtendsIn, SpinorFillTypes{SpinorFillType::one}, true, typeOfComparisionIn) {};
	EvenOddSpinorTestParameters() : SpinorTestParameters() {};
};

//todo: need children for evenOdd and nonEvenOdd. Then, one can have a common "size" member instead of spinorfieldsize and evenOddSpinorfieldsize. This should simplify the usage a lot!
class SpinorTester : public KernelTester {
public: //@todo: make these protected!
	SpinorTester(std::string kernelName, const ParameterCollection,	const SpinorTestParameters & ); //@todo: this must go away in the end!
	SpinorTester(std::string kernelName, const ParameterCollection,	const SpinorTestParameters &, const bool, const size_t );
	SpinorTester(std::string kernelName, const ParameterCollection,	const SpinorTestParameters &, const bool, const size_t, const ReferenceValues );
	
protected:
	spinor * createSpinorfield( SpinorFillType ); // @todo: this always create a nonEo field!!!
	spinor * createSpinorfieldWithOnesAndZerosDependingOnSiteParity(const bool fillEvenSites);
	spinor * createSpinorfieldWithOnesAndMinusOneForGamma5Use(size_t numberOfElements);	
	//todo: these should not be visible here, but be accessible via a fillType
	void fillTwoSpinorBuffers(const hardware::buffers::Spinor * in1, const hardware::buffers::Spinor * in2, int seed = 123456); //this is used in the molecular dynamics test
	void fillTwoSpinorBuffersDependingOnParity(const hardware::buffers::Spinor * in1, const hardware::buffers::Spinor * in2);
	void fillTwoSpinorfieldsDependingOnParity(spinor * sf_in1, spinor * sf_in2, int size);
	void fillTwoSpinorfieldsWithRandomNumbers(spinor * sf_in1, spinor * sf_in2, int size, int seed = 123456);

	hmc_float count_sf(spinor * in, int size);
	hmc_float calc_var_sf(spinor * in, int size, hmc_float sum);
	void calcSquarenormAndStoreAsKernelResult(const hardware::buffers::Plain<spinor> * in);
	void calcSquarenormEvenOddAndStoreAsKernelResult(const hardware::buffers::Spinor * in);
	
	int getSpinorfieldSize(const LatticeExtents latticeExtendsIn) const { return calculateSpinorfieldSize(latticeExtendsIn); } ;
	int getEvenOddSpinorfieldSize(const LatticeExtents latticeExtendsIn) const { return calculateEvenOddSpinorfieldSize(latticeExtendsIn); } ;

	const hardware::code::Spinors * code;

	hardware::buffers::Plain<double> * doubleBuffer;
	
	//todo: most of these must go away!
	size_t spinorfieldElements;
	size_t spinorfieldEvenOddElements;
	bool evenOrOdd;
	const bool isEvenOdd;
	const size_t elements;
};

struct NonEvenOddSpinorTester : public SpinorTester
{
	NonEvenOddSpinorTester(const std::string kernelName, const ParameterCollection pC, const SpinorTestParameters & tP) :
		SpinorTester(kernelName, pC, tP, false, getSpinorfieldSize(tP.latticeExtents)) {};
	NonEvenOddSpinorTester(const std::string kernelName, const ParameterCollection pC, const SpinorTestParameters & tP, const ReferenceValues (*rV) (const int, const SpinorFillTypes sF)) :
		SpinorTester(kernelName, pC, tP, false, getSpinorfieldSize(tP.latticeExtents), rV(getSpinorfieldSize(tP.latticeExtents), tP.fillTypes)) {};
	NonEvenOddSpinorTester(const std::string kernelName, const ParameterCollection pC, const SpinorTestParameters & tP, const ReferenceValues & rV) :
		SpinorTester(kernelName, pC, tP, false, getSpinorfieldSize(tP.latticeExtents), rV) {};
};

struct EvenOddSpinorTester : public SpinorTester
{
	EvenOddSpinorTester(const std::string kernelName, const ParameterCollection pC, const SpinorTestParameters & tP) :
		SpinorTester(kernelName, pC, tP, true, getEvenOddSpinorfieldSize(tP.latticeExtents)) {};
	EvenOddSpinorTester(const std::string kernelName, const ParameterCollection pC, const SpinorTestParameters & tP, const ReferenceValues & rV) :
		SpinorTester(kernelName, pC, tP, true, getEvenOddSpinorfieldSize(tP.latticeExtents), rV) {};
};

