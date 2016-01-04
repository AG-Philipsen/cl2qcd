/** @file
 * Testcases for the hardware::buffers::Buffer class
 *
 * Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
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

#include "buffer.hpp"
#include "plain.hpp"
#include "../system.hpp"
#include "../device.hpp"
#include "../../meta/type_ops.hpp"
#include <stdexcept>
#include "../interfaceMockups.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hardware::buffers::Buffer
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace hardware;
	using namespace hardware::buffers;

	const hardware::HardwareParametersMockup hardwareParameters(4,4);
	const hardware::code::OpenClKernelParametersMockup kernelParameters(4,4);
	hardware::System system( hardwareParameters, kernelParameters );
	const std::vector<Device*>& devices = system.get_devices();
	for(Device * device : devices)
	{
		Buffer dummy(sizeof(float), device);
		BOOST_REQUIRE_EQUAL(dummy.get_bytes(), sizeof(float));
		const cl_mem * tmp = dummy;
		BOOST_REQUIRE(tmp);
		BOOST_REQUIRE(*tmp);
		BOOST_REQUIRE_EQUAL(device->get_id(), dummy.get_device()->get_id());
	}
}

BOOST_AUTO_TEST_CASE(copy)
{
	// normal buffer doesn't have input/output -> only test exception if different size

	using namespace hardware;
	using namespace hardware::buffers;

	const hardware::HardwareParametersMockup hardwareParameters(4,4);
	const hardware::code::OpenClKernelParametersMockup kernelParameters(4,4);
	hardware::System system( hardwareParameters, kernelParameters );
	const std::vector<Device*>& devices = system.get_devices();
	for(Device * device : devices)
	{
		Buffer dummy(sizeof(float), device);
		BOOST_REQUIRE_EQUAL(dummy.get_bytes(), sizeof(float));
		Buffer dummy2(2 * dummy.get_bytes(), device);
		Buffer dummy3(2 * dummy.get_bytes(), device);

		BOOST_CHECK_THROW(copyData(&dummy, &dummy2), std::invalid_argument);
		copyData(&dummy2, &dummy3);
		BOOST_CHECK_THROW(copyData(&dummy3, &dummy), std::invalid_argument);
	}
}

BOOST_AUTO_TEST_CASE(clear)
{
	// normal buffer doesn't have input/output -> only test exception if different size

	using namespace hardware;
	using namespace hardware::buffers;

	const hardware::HardwareParametersMockup hardwareParameters(4,4);
	const hardware::code::OpenClKernelParametersMockup kernelParameters(4,4);
	hardware::System system( hardwareParameters, kernelParameters );
	const std::vector<Device*>& devices = system.get_devices();
	for(Device * device : devices)
	{
		{
			const size_t SIZE = 1024;
			float host[SIZE];
			fill(host, SIZE);
			Plain<float> dummy(SIZE, device);
			dummy.load(host);
			dummy.clear();
			dummy.dump(host);
			for(size_t i = 0; i < SIZE; ++i) {
				BOOST_REQUIRE_EQUAL(host[i], 0.f);
			}
		}

		{
			const size_t SIZE = 71;
			char host[SIZE];
			fill(host, SIZE);
			Plain<char> dummy(SIZE, device);
			dummy.load(host);
			dummy.clear();
			dummy.dump(host);
			for(size_t i = 0; i < SIZE; ++i) {
				BOOST_REQUIRE_EQUAL(host[i], 0);
			}
		}
	}
}
