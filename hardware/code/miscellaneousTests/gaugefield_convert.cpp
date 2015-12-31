/*
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

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE gaugefield_convert
#include <boost/test/unit_test.hpp>

#include "../../../meta/type_ops.hpp" //@todo: move the Matrixsu3 fcts. from here to a different place

#include "../mockups.hpp"
#include "../GaugefieldTester.hpp"

void test(const hardware::System& system, const int seed, const LatticeExtents lE)
{
	const size_t NUM_ELEMENTS = calculateGaugefieldSize(lE);
	for(auto device: system.get_devices())
	{
		Matrixsu3 * const in = new Matrixsu3[NUM_ELEMENTS];
		Matrixsu3 * const out = new Matrixsu3[NUM_ELEMENTS];
		fill(in, NUM_ELEMENTS, seed);
		fill(out, NUM_ELEMENTS, seed + NUM_ELEMENTS);
		hardware::buffers::SU3 buffer(NUM_ELEMENTS, device);

		auto code = device->getGaugefieldCode();
		code->importGaugefield(&buffer, in);
		code->exportGaugefield(out, &buffer);

		BOOST_CHECK_EQUAL_COLLECTIONS(in, in + NUM_ELEMENTS, out, out + NUM_ELEMENTS);

		delete[] in;
		delete[] out;
	}
}

BOOST_AUTO_TEST_CASE(CPU)
{
	LatticeExtents lE{4,4};
	hardware::HardwareParametersMockup hardwareParameters(lE.ns,lE.nt);
	hardware::code::OpenClKernelParametersMockup kernelParameters(lE.ns,lE.nt);
	hardware::System system(hardwareParameters, kernelParameters);

	test(system, 1, lE);
	test(system, 14, lE);
	test(system, 21, lE);
}

BOOST_AUTO_TEST_CASE(GPU)
{
	BOOST_ERROR("not implemented");
	//@todo: This is in principle the same as for the CPU, one has to see how one switches between GPU and CPU in general first...
}



