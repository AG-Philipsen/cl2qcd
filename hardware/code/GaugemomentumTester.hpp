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

#ifndef GAUGEMOMENTUM_TESTER_HPP_
#define GAUGEMOMENTUM_TESTER_HPP_

#include "kernelTester.hpp"
#include "gaugemomentum.hpp"
#include "../../host_functionality/host_random.h"
#include "prng.hpp"
#include "SpinorTester.hpp"

int calculateGaugemomentumSize(const LatticeExtents latticeExtentsIn) noexcept;
int calculateAlgebraSize(const LatticeExtents latticeExtentsIn) noexcept;

struct GaugemomentumTestParameters: public TestParameters
{
	GaugemomentumTestParameters(const LatticeExtents latticeExtendsIn) :
		TestParameters(latticeExtendsIn), fillType(SpinorFillType::one)
	{};

	const SpinorFillType fillType;
};

class GaugemomentumTester : public KernelTester
{
public:
	GaugemomentumTester(const std::string kernelName, const ParameterCollection pC, const ReferenceValues rV, const GaugemomentumTestParameters tP);
	virtual ~GaugemomentumTester();
protected:
	const hardware::code::Gaugemomentum * code;
	hardware::buffers::Plain<double> * doubleBuffer;
	void calcSquarenormAndStoreAsKernelResult(const hardware::buffers::Gaugemomentum * in, int index = 0);
};


struct GaugemomentumCreator
{
	GaugemomentumCreator(size_t numberOfElementsIn): numberOfElements(numberOfElementsIn){};
	enum Filltype {one, zero};
	double * createGaugemomentumBasedOnFilltype(const Filltype filltype = one);
	void fill_with_one(double * sf_in);
	void fill_with_zero(double * sf_in);
	double count_gm(ae * ae_in, int size);
	double calc_var(double in, double mean);
	double calc_var_gm(ae * ae_in, int size, double sum);

	size_t numberOfElements;

};

#endif
