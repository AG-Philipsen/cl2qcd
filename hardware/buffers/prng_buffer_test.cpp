/** @file
 * Testcases for the hardware::buffers::Buffer class
 *
 * Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>,
 * 		2015 Christopher Pinke
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

#include "prng_buffer.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hardware::buffers::PRNGBuffer
#include <boost/test/unit_test.hpp>

#include "../system.hpp"
#include "../interfaceMockups.hpp"

BOOST_AUTO_TEST_CASE(get_prng_buffer_size)
{
	using namespace hardware;
	using namespace hardware::buffers;

	const hardware::HardwareParametersMockup hardwareParameters(4,4);
	const hardware::code::OpenClKernelParametersMockup kernelParameters(4,4);
	hardware::System system( hardwareParameters, kernelParameters );
	for(Device * device : system.get_devices())
	{
		int elems = hardware::buffers::get_prng_buffer_size(device, hardwareParameters.useSameRandomNumbers());

		BOOST_CHECK_GT(elems, 0);
		BOOST_CHECK_LT(elems, 1e6);
	}
}

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace hardware;
	using namespace hardware::buffers;

	const hardware::HardwareParametersMockup hardwareParameters(4,4);
	const hardware::code::OpenClKernelParametersMockup kernelParameters(4,4);
	hardware::System system( hardwareParameters, kernelParameters );
	for(Device * device : system.get_devices())
	{
		PRNGBuffer dummy(device, hardwareParameters.useSameRandomNumbers());
		BOOST_CHECK_EQUAL(dummy.get_elements(), hardware::buffers::get_prng_buffer_size(device, hardwareParameters.useSameRandomNumbers()));
		const cl_mem * tmp = dummy;
		BOOST_CHECK(tmp);
		BOOST_CHECK(*tmp);
	}
}

BOOST_AUTO_TEST_CASE(import_export)
{
	using namespace hardware;
	using namespace hardware::buffers;

	const hardware::HardwareParametersMockup hardwareParameters(4,4);
	const hardware::code::OpenClKernelParametersMockup kernelParameters(4,4);
	hardware::System system( hardwareParameters, kernelParameters );
	for(Device * device : system.get_devices())
	{
		PRNGBuffer buffer(device, hardwareParameters.useSameRandomNumbers());
		int elems = buffer.get_elements();
		PRNGBuffer::prng_state_t * in = new PRNGBuffer::prng_state_t[elems];
		// TODO fill with random data
		buffer.load(in);
		PRNGBuffer::prng_state_t * out = new PRNGBuffer::prng_state_t[elems];
		buffer.dump(out);

		BOOST_CHECK_EQUAL_COLLECTIONS(reinterpret_cast<uint64_t*>(in), reinterpret_cast<uint64_t*>(in + elems), reinterpret_cast<uint64_t*>(out), reinterpret_cast<uint64_t*>(out + elems));

		delete[] out;
		delete[] in;
	}
}
