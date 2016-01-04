/** @file
 * Testcases for the hardware::buffers::Matrix3x3 class
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

#include "3x3.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hardware::buffers::Matrix3x3
#include <boost/test/unit_test.hpp>

#include "../system.hpp"
#include "../interfaceMockups.hpp"

BOOST_AUTO_TEST_CASE(initialization)
{
	using namespace hardware;

	const hardware::HardwareParametersMockup hardwareParameters(4,4);
	const hardware::code::OpenClKernelParametersMockup kernelParameters(4,4);
	hardware::System system( hardwareParameters, kernelParameters );
	for(Device * device : system.get_devices())
	{
		hardware::buffers::Matrix3x3 dummy(system.getHardwareParameters()->getLatticeVolume() * NDIM, device);
		const cl_mem * tmp = dummy;
		BOOST_CHECK(tmp);
		BOOST_CHECK(*tmp);
	}
}
