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
typedef std::vector<hmc_float> RealNumbers;
typedef size_t NumberOfSpinors;

int calculateSpinorfieldSize(LatticeExtents latticeExtendsIn) noexcept;
int calculateEvenOddSpinorfieldSize(const LatticeExtents latticeExtendsIn) noexcept;

hmc_float count_sf(spinor * in, int size);
hmc_float calc_var(hmc_float in, hmc_float mean);
hmc_float calc_var_sf(spinor * in, int size, hmc_float sum);

struct SpinorTestParameters: public virtual TestParameters
{
	SpinorTestParameters(const LatticeExtents latticeExtendsIn) :
		TestParameters(latticeExtendsIn), fillTypes(SpinorFillType::one) {};
	SpinorTestParameters(const LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn) :
		TestParameters(latticeExtendsIn), fillTypes(fillTypesIn) {};

	const SpinorFillTypes fillTypes;
};

class SpinorTester: public KernelTester
{
public:
	SpinorTester(std::string kernelName, const ParameterCollection, const SpinorTestParameters &, const ReferenceValues );
protected:
	void calcSquarenormAndStoreAsKernelResult(const hardware::buffers::Plain<spinor> * in);
	void calcSquarenormEvenOddAndStoreAsKernelResult(const hardware::buffers::Spinor * in);

	const hardware::code::Spinors * code;
	hardware::buffers::Plain<double> * doubleBuffer;
};

struct SpinorfieldCreator
{
	SpinorfieldCreator(const size_t numberOfElementsIn): numberOfElements(numberOfElementsIn) {};
	spinor * createSpinorfield(SpinorFillType);

	size_t numberOfElements;
};

struct NonEvenOddSpinorfieldCreator : public SpinorfieldCreator
{
	NonEvenOddSpinorfieldCreator(LatticeExtents lE): SpinorfieldCreator(calculateSpinorfieldSize(lE)), latticeExtents(lE){};
	spinor * createSpinorfieldWithOnesAndZerosDependingOnSiteParity(const bool fillEvenSites);
	LatticeExtents latticeExtents;
};

struct EvenOddSpinorfieldCreator : public SpinorfieldCreator
{
	EvenOddSpinorfieldCreator(LatticeExtents lE): SpinorfieldCreator(calculateEvenOddSpinorfieldSize(lE)), latticeExtents(lE){};
	//todo: these should not be visible here, but be accessible via a fillType
	void fillTwoSpinorBuffers(const hardware::buffers::Spinor * in1, const SpinorFillType fillTypeIn1, const hardware::buffers::Spinor * in2, const SpinorFillType fillTypeIn2);
	void fillTwoSpinorfieldsDependingOnParity(spinor * sf_in1, spinor * sf_in2, int size);
	void fillTwoSpinorBuffersDependingOnParity(const hardware::buffers::Spinor * in1, const hardware::buffers::Spinor * in2);
	LatticeExtents latticeExtents;
};

struct NonEvenOddSpinorTester : public SpinorTester
{
	NonEvenOddSpinorTester(const std::string kernelName, const ParameterCollection pC, const SpinorTestParameters & tP, const ReferenceValues & rV) :
		SpinorTester(kernelName, pC, tP, rV) {};
};

struct EvenOddSpinorTester : public SpinorTester
{
	EvenOddSpinorTester(const std::string kernelName, const ParameterCollection pC, const SpinorTestParameters & tP, const ReferenceValues & rV) :
		SpinorTester(kernelName, pC, tP, rV) {};
};

