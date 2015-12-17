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

//@todo: move to common place
//@todo: rename to FillType
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

struct SpinorTestParameters: public virtual TestParameters
{
	SpinorTestParameters(const LatticeExtents latticeExtendsIn) :
		TestParameters(latticeExtendsIn), fillTypes(SpinorFillType::one) {};
	SpinorTestParameters(const LatticeExtents latticeExtendsIn, const ComparisonType typeOfComparisionIn) :
		TestParameters(latticeExtendsIn, typeOfComparisionIn), fillTypes(SpinorFillType::one) {};
	SpinorTestParameters(const LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn) :
		TestParameters(latticeExtendsIn), fillTypes(fillTypesIn) {};

	const SpinorFillTypes fillTypes;
};

struct SpinorTester : public KernelTester
{
	SpinorTester(std::string kernelName, const ParameterCollection,	const SpinorTestParameters &, const size_t, const ReferenceValues );
protected:
	spinor * createSpinorfield( SpinorFillType );
	//@todo: check if all of these fcts. are actually used
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
	const size_t elements;
};

struct NonEvenOddSpinorTester : public SpinorTester
{
	NonEvenOddSpinorTester(const std::string kernelName, const ParameterCollection pC, const SpinorTestParameters & tP, const ReferenceValues & rV) :
		SpinorTester(kernelName, pC, tP, getSpinorfieldSize(tP.latticeExtents), rV) {};
};

struct EvenOddSpinorTester : public SpinorTester
{
	EvenOddSpinorTester(const std::string kernelName, const ParameterCollection pC, const SpinorTestParameters & tP, const ReferenceValues & rV) :
		SpinorTester(kernelName, pC, tP, getEvenOddSpinorfieldSize(tP.latticeExtents), rV) {};
};

