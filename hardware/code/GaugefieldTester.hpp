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

#include "kernelTester.hpp"

enum GaugefieldFillType {cold = 1, nonTrivial};

struct GaugefieldTestParameters : public TestParameters
{
	GaugefieldFillType fillType;

	GaugefieldTestParameters(std::vector<double> referenceValueIn, const LatticeExtents latticeExtentsIn, GaugefieldFillType fillTypeIn):
		TestParameters(referenceValueIn, latticeExtentsIn), fillType( fillTypeIn ) {}
};

class GaugefieldTester : public KernelTester {
public:
	GaugefieldTester(std::string, const ParameterCollection &, const GaugefieldTestParameters);
	~GaugefieldTester();

protected:
	const int numberOfElements;
	const hardware::code::Gaugefield * code;
	const hardware::buffers::SU3 * gaugefieldBuffer;
};





