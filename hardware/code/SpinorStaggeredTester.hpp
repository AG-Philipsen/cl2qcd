/** @file
 * SpinorStaggeredTester class
 *
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
#ifndef _HARDWARE_CODE_SPINOR_STAGGERED_TESTER_
#define _HARDWARE_CODE_SPINOR_STAGGERED_TESTER_

#include "kernelTester.hpp"
#include "spinors_staggered.hpp"
#include "SpinorTester.hpp"

//Utilities methods
std::vector<hmc_float> reals_from_su3vec(su3vec v);
hmc_float count_sf(su3vec * in, int size);
hmc_float squareNorm(su3vec * in, int size);
hmc_float calc_var_sf(su3vec * in, int size, hmc_float sum);

struct SpinorStaggeredTestParameters: public virtual TestParameters
{
	SpinorStaggeredTestParameters(const LatticeExtents latticeExtendsIn) :
		TestParameters(latticeExtendsIn), fillTypes(SpinorFillType::one) {};
	SpinorStaggeredTestParameters(const LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn) :
		TestParameters(latticeExtendsIn), fillTypes(fillTypesIn) {};

	const SpinorFillTypes fillTypes;
};

class SpinorStaggeredTester : public KernelTester
{
public:
	SpinorStaggeredTester(const std::string kernelName, const ParameterCollection, const SpinorStaggeredTestParameters &, const ReferenceValues);
protected:
	void calcSquarenormAndStoreAsKernelResult(const hardware::buffers::Plain<su3vec> * in);
	void calcSquarenormEvenOddAndStoreAsKernelResult(const hardware::buffers::SU3vec * in);

	const hardware::code::Spinors_staggered * code;
	hardware::buffers::Plain<double> * doubleBuffer;

};

struct SpinorStaggeredfieldCreator
{
	SpinorStaggeredfieldCreator(const size_t numberOfElementsIn): numberOfElements(numberOfElementsIn) {};
	su3vec * createSpinorfield(SpinorFillType);

	size_t numberOfElements;
};

struct NonEvenOddSpinorStaggeredfieldCreator : public SpinorStaggeredfieldCreator
{
	NonEvenOddSpinorStaggeredfieldCreator(const LatticeExtents lE): SpinorStaggeredfieldCreator(calculateSpinorfieldSize(lE)), latticeExtents(lE){};
	su3vec * createSpinorfieldWithOnesAndZerosDependingOnSiteParity(const bool fillEvenSites);
	LatticeExtents latticeExtents;
};

struct EvenOddSpinorStaggeredfieldCreator : public SpinorStaggeredfieldCreator
{
	EvenOddSpinorStaggeredfieldCreator(const LatticeExtents lE): SpinorStaggeredfieldCreator(calculateEvenOddSpinorfieldSize(lE)), latticeExtents(lE){};
	su3vec * createSpinorfieldEvenOddWithOnesAndZerosDependingOnSiteParity(const bool fillEvenSites);
	LatticeExtents latticeExtents;
};

struct NonEvenOddSpinorStaggeredTester : public SpinorStaggeredTester
{
	NonEvenOddSpinorStaggeredTester(const std::string kernelName, const ParameterCollection pC, const SpinorStaggeredTestParameters & tP, const ReferenceValues & rV) :
		SpinorStaggeredTester(kernelName, pC, tP, rV) {};
};

struct EvenOddSpinorStaggeredTester : public SpinorStaggeredTester
{
	EvenOddSpinorStaggeredTester(const std::string kernelName, const ParameterCollection pC, const SpinorStaggeredTestParameters & tP, const ReferenceValues & rV) :
		SpinorStaggeredTester(kernelName, pC, tP, rV) {};
};

#endif // _HARDWARE_CODE_SPINOR_STAGGERED_TESTER_
