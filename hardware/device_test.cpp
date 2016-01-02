/** @file
 * Testcases for the hardware::Device class
 *
 * Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>,
 * 	2015 Christopher Pinke <pinke@th.physik.uni-frankfurt.de>
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
#define BOOST_TEST_MODULE hardware::Device
#include <boost/test/unit_test.hpp>

#include "system.hpp"
#include "device.hpp"
#include "hardwareTestUtilities.hpp"
#include "interfaceMockups.hpp"

void querrySomeInformationsFromDevice( const hardware::Device * device )
{
	BOOST_REQUIRE_NE(device->get_preferred_local_thread_num(), 0);
	BOOST_REQUIRE_NE(device->get_preferred_global_thread_num(), 0);
	const size_t recommended_stride = device->recommendStride(1024, 16, 2);
	BOOST_REQUIRE_GE(recommended_stride, 1024);
}

BOOST_AUTO_TEST_CASE(initialization)
{
	const hardware::HardwareParametersMockup hardwareParameters(4,4);
	const hardware::code::OpenClKernelParametersMockup kernelParameters(4,4);
	hardware::System system( hardwareParameters, kernelParameters );
	atLeastOneDeviceMustExistForSanityOfSystem( &system );

	for(const hardware::Device * device : system.get_devices())
	{
		querrySomeInformationsFromDevice(device);
	}
}

BOOST_AUTO_TEST_CASE(compile)
{
	const hardware::HardwareParametersMockup hardwareParameters(4,4);
	const hardware::code::OpenClKernelParametersMockup kernelParameters(4,4);
	hardware::System system( hardwareParameters, kernelParameters );
	atLeastOneDeviceMustExistForSanityOfSystem( &system );

	for(const hardware::Device * device : system.get_devices())
	{
		device->createKernel("foo", "") << "../hardware/device_test.cl";
	}
}
