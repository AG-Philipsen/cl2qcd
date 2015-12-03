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

struct SpinorStaggeredTestParameters: public virtual TestParameters
{
	SpinorStaggeredTestParameters(const LatticeExtents latticeExtendsIn) :
		TestParameters(latticeExtendsIn), fillTypes(SpinorFillType::one) {};
	SpinorStaggeredTestParameters(const LatticeExtents latticeExtendsIn, const int typeOfComparisionIn) :
		TestParameters(latticeExtendsIn, typeOfComparisionIn), fillTypes(SpinorFillType::one) {};
	SpinorStaggeredTestParameters(const LatticeExtents latticeExtendsIn, const SpinorFillTypes fillTypesIn) :
		TestParameters(latticeExtendsIn), fillTypes(fillTypesIn) {};

	const SpinorFillTypes fillTypes;
};

class SpinorStaggeredTester : public KernelTester {
public:
	SpinorStaggeredTester(const std::string kernelName, const ParameterCollection, const SpinorStaggeredTestParameters &, const size_t, const ReferenceValues);
protected:
	//Methods (protected for inheritance resons)
	su3vec * createSpinorfield(SpinorFillType);
	su3vec * createSpinorfieldWithOnesAndZerosDependingOnSiteParity(const bool fillEvenSites);
	su3vec * createSpinorfieldEvenOddWithOnesAndZerosDependingOnSiteParity(const bool fillEvenSites);

	void fill_with_one_eo(su3vec * sf_in, int size, bool eo);
	hmc_float count_sf_eo(su3vec * sf_in, int size, bool eo);
	
	//These methods are used to produce files for the Reference Code (D'Elia et al)
	void print_staggeredfield_to_textfile(std::string outputfile, su3vec * sf);
	void print_staggeredfield_eo_to_textfile(std::string outputfile, su3vec * sf);
	//Utilities methods
	std::vector<hmc_float> reals_from_su3vec(su3vec v);

	hmc_float count_sf(su3vec * in, int size);
	hmc_float calc_var(hmc_float in, hmc_float mean);
	hmc_float calc_var_sf(su3vec * in, int size, hmc_float sum);
	void calcSquarenormAndStoreAsKernelResult(const hardware::buffers::Plain<su3vec> * in);
	void calcSquarenormEvenOddAndStoreAsKernelResult(const hardware::buffers::SU3vec * in);

	int getSpinorfieldSize(const LatticeExtents latticeExtendsIn) const { return calculateSpinorfieldSize(latticeExtendsIn); } ;
	int getEvenOddSpinorfieldSize(const LatticeExtents latticeExtendsIn) const { return calculateEvenOddSpinorfieldSize(latticeExtendsIn); } ;

	const hardware::code::Spinors_staggered * code;
	hardware::buffers::Plain<double> * doubleBuffer;
	const size_t elements;
};

struct NonEvenOddSpinorStaggeredTester : public SpinorStaggeredTester
{
	NonEvenOddSpinorStaggeredTester(const std::string kernelName, const ParameterCollection pC, const SpinorStaggeredTestParameters & tP, const ReferenceValues & rV) :
		SpinorStaggeredTester(kernelName, pC, tP, getSpinorfieldSize(tP.latticeExtents), rV) {};
};

struct EvenOddSpinorStaggeredTester : public SpinorStaggeredTester
{
	EvenOddSpinorStaggeredTester(const std::string kernelName, const ParameterCollection pC, const SpinorStaggeredTestParameters & tP, const ReferenceValues & rV) :
		SpinorStaggeredTester(kernelName, pC, tP, getEvenOddSpinorfieldSize(tP.latticeExtents), rV) {};
};

#endif // _HARDWARE_CODE_SPINOR_STAGGERED_TESTER_
