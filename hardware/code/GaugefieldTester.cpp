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

#include "GaugefieldTester.hpp"
//#include "../../host_functionality/host_operations_gaugefield.h"

Matrixsu3 unit_matrixsu3()
{
	Matrixsu3 out;
	out.e00.re = 1.;
	out.e00.im = 0.;
	out.e01.re = 0.;
	out.e01.im = 0.;
	out.e02.re = 0.;
	out.e02.im = 0.;

	out.e10.re = 0.;
	out.e10.im = 0.;
	out.e11.re = 1.;
	out.e11.im = 0.;
	out.e12.re = 0.;
	out.e12.im = 0.;

	out.e20.re = 0.;
	out.e20.im = 0.;
	out.e21.re = 0.;
	out.e21.im = 0.;
	out.e22.re = 1.;
	out.e22.im = 0.;

	return out;
}

Matrixsu3 nonTrivialSu3Matrix()
{
	Matrixsu3 out;
	out.e00.re = .130189;
	out.e00.im = .260378;
	out.e01.re = .260378;
	out.e01.im = .390567;
	out.e02.re = .520756;
	out.e02.im = .650945;

	out.e10.re = .572742;
	out.e10.im = .403041;
	out.e11.re = .371222;
	out.e11.im = .321726;
	out.e12.re = -.449002;
	out.e12.im = -.258088;

	out.e20.re = 1.11022e-16;
	out.e20.im = .651751;
	out.e21.re = .0271563;
	out.e21.im = -.733219;
	out.e22.re = -.0271563;
	out.e22.im = .190094;

	return out;
}

GaugefieldTester::GaugefieldTester(std::string kernelName, const ParameterCollection & parameterCollection, const GaugefieldTestParameters testParams, const ReferenceValues rV):
	KernelTester(kernelName, parameterCollection.hardwareParameters, parameterCollection.kernelParameters, testParams, rV)
{
	code = device->getGaugefieldCode();
}

int calculateGaugefieldSize(const LatticeExtents latticeExtentsIn) noexcept
{
	return 	latticeExtentsIn.getLatticeVolume() * NDIM;
}

const Matrixsu3* GaugefieldCreator::createGaugefield(const GaugefieldFillType fillTypeIn)
{
	Matrixsu3 * tmp = new Matrixsu3[numberOfElements];
	for(size_t i = 0; i < (size_t)numberOfElements; ++i)
	{
		switch (fillTypeIn)
		{
		case GaugefieldFillType::cold:
			tmp[i] = unit_matrixsu3();
			break;
		case GaugefieldFillType::nonTrivial:
			tmp[i] = nonTrivialSu3Matrix();
			break;
		default:
			BOOST_ERROR("No valid GaugefieldFillType specified");
		}
	}
	return tmp;
}
