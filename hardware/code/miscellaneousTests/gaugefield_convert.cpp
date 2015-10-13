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

#include "../../../meta/util.hpp"
#include "../../../meta/type_ops.hpp"
#include "../../device.hpp"
#include "../../system.hpp"
#include "../gaugefield.hpp"

void test(const hardware::System& system, const int seed)
{
	const auto & params = system.get_inputparameters();
	const size_t NUM_ELEMENTS = meta::get_vol4d(params) * NDIM;
for(auto device: system.get_devices()) {
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
	const char* _params[] = {"foo", "--use_gpu=false"};
	meta::Inputparameters params(2, _params);
	hardware::System system(params);

	test(system, 1);
	test(system, 14);
	test(system, 21);
}

BOOST_AUTO_TEST_CASE(GPU)
{
	const char* _params[] = {"foo", "--use_cpu=false"};
	meta::Inputparameters params(2, _params);
	hardware::System system(params);

	test(system, 1);
	test(system, 14);
	test(system, 21);
}
