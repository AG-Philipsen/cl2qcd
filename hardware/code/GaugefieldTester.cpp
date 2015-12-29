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

#include "../../host_functionality/host_operations_gaugefield.h"

void setGaugefield(Matrixsu3 * field, size_t elems, const GaugefieldFillType fillTypeIn)
{
	for(size_t i = 0; i < elems; ++i)
	{
		switch (fillTypeIn)
		{
		case GaugefieldFillType::cold:
			field[i] = unit_matrixsu3();
			break;
		case GaugefieldFillType::nonTrivial:
			field[i] = nonTrivialSu3Matrix();
			break;
		default:
			BOOST_ERROR("No valid GaugefieldFillType specified");
		}
	}
}

const Matrixsu3* createGaugefield(const int numberOfElements, const GaugefieldFillType fillTypeIn)
{
	Matrixsu3 * tmp = new Matrixsu3[numberOfElements];
	setGaugefield( tmp, numberOfElements, fillTypeIn);
	return tmp;
}

GaugefieldTester::GaugefieldTester(std::string kernelName, const ParameterCollection & parameterCollection, const GaugefieldTestParameters testParams, const ReferenceValues rV):
	KernelTester(kernelName, parameterCollection.hardwareParameters, parameterCollection.kernelParameters, testParams, rV),
	numberOfElements(parameterCollection.hardwareParameters.getLatticeVolume() * NDIM)
{
	gaugefieldBuffer = new hardware::buffers::SU3( numberOfElements, device);
	const Matrixsu3 * gf_host = createGaugefield(numberOfElements, testParams.fillType);
	device->getGaugefieldCode()->importGaugefield(gaugefieldBuffer, gf_host);
	delete[] gf_host;

	code = device->getGaugefieldCode();
}

GaugefieldTester::~GaugefieldTester()
{
	delete gaugefieldBuffer;
}



