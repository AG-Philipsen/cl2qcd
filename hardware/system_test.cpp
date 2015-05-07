/** @file
 * Testcases for the hardware::System class
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hardware::System
#include <boost/test/unit_test.hpp>

#include "system.hpp"
#include "device.hpp"

const char * dummyRuntimeArguments[] = {"foo"};
const int dummyNumberOfRuntimeArguments = 1;

void broadcastMessage_warn(const std::string message)
{
	logger.warn() << message;
	BOOST_TEST_MESSAGE( message );
}

void broadcastMessage_fatal(const std::string message)
{
	logger.fatal() << message;
	BOOST_TEST_MESSAGE( message );
}

void failTest()
{
	BOOST_REQUIRE_EQUAL(true, false);
}

BOOST_AUTO_TEST_CASE(initialization)
{
	meta::Inputparameters params(dummyNumberOfRuntimeArguments, dummyRuntimeArguments);
	hardware::System system( params );
	BOOST_REQUIRE_EQUAL(&system.get_inputparameters(), &params);
}

BOOST_AUTO_TEST_SUITE(systemSanity)

	void atLeastOneDeviceMustExistForSanityOfSystem(const hardware::System * system)
	{
		BOOST_REQUIRE_GE(system->get_devices().size(), 1);
	}

	void allDevicesMustSupportDoublePrecisionForSanityOfSystem(const hardware::System * system)
	{
		for(hardware::Device * device : system->get_devices())
		{
			BOOST_REQUIRE_EQUAL(device->is_double_supported(), true);
		}
	}

	BOOST_AUTO_TEST_CASE(staticCastToPlatform)
	{
		meta::Inputparameters params(dummyNumberOfRuntimeArguments, dummyRuntimeArguments);
		hardware::System system(params);
		BOOST_REQUIRE(static_cast<const cl_context>(system));
	}

	BOOST_AUTO_TEST_CASE(enoughDevicesExist)
	{
		meta::Inputparameters params(dummyNumberOfRuntimeArguments, dummyRuntimeArguments);
		hardware::System system(params);
		atLeastOneDeviceMustExistForSanityOfSystem( &system );
	}

	BOOST_AUTO_TEST_CASE(allDevicesHaveDoubleSupport)
	{
		meta::Inputparameters params(dummyNumberOfRuntimeArguments, dummyRuntimeArguments);
		hardware::System system(params);
		allDevicesMustSupportDoublePrecisionForSanityOfSystem( &system );
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(devices)

	void checkThatOnlySpecifiedNumberOfDevicesIsInitialized(const int totalNumberOfDevicesInSystem)
	{
		for( int desiredNumberOfDevices = 1; desiredNumberOfDevices <= totalNumberOfDevicesInSystem; desiredNumberOfDevices ++)
		{
			const std::string tmp = ("--num_dev=" + std::to_string( desiredNumberOfDevices ) );
			const char * _params[] = {"foo", tmp.c_str()};
			meta::Inputparameters params(2, _params);
			hardware::System system(params);
			BOOST_REQUIRE_EQUAL( system.get_devices().size() , desiredNumberOfDevices);
		}
	}

	BOOST_AUTO_TEST_CASE(setNumberOfDevicesByCommandLine)
	{
		meta::Inputparameters params(dummyNumberOfRuntimeArguments, dummyRuntimeArguments);
		hardware::System system(params);
		checkThatOnlySpecifiedNumberOfDevicesIsInitialized( system.get_devices().size());
	}

	void checkThatOnlySpecifiedDeviceIsInitialized(const int totalNumberOfDevicesInSystem)
	{
		for( int desiredDevice = 0; desiredDevice < totalNumberOfDevicesInSystem; desiredDevice ++)
		{
			const std::string tmp = ("--device=" + std::to_string( desiredDevice ) );
			const char * _params[] = {"foo", tmp.c_str()};
			meta::Inputparameters params(2, _params);
			try {
				hardware::System system(params);
				BOOST_REQUIRE_EQUAL( system.get_devices().size(), 1);
			} catch(std::invalid_argument) {
				// device might not support double-precision
			}
		}
	}

	BOOST_AUTO_TEST_CASE(setDeviceByCommandLine)
	{
		meta::Inputparameters params(dummyNumberOfRuntimeArguments, dummyRuntimeArguments);
		hardware::System system(params);
		checkThatOnlySpecifiedDeviceIsInitialized( system.get_devices().size());
	}

	bool checkIfNoOpenCLDevicesWereFound( const hardware::OpenclException exception)
	{
		return exception.errorCode == -1;
	}

	void disableSpecificDeviceTypeByCommandLine( const char * _parameters[], const cl_device_type device_type)
	{
		meta::Inputparameters params(2, _parameters);
		try {
			hardware::System system(params);
			for(hardware::Device * device : system.get_devices())
			{
				BOOST_REQUIRE_NE(device->get_device_type(), device_type);
			}
		}
		catch(hardware::OpenclException exception)
		{
			if ( checkIfNoOpenCLDevicesWereFound( exception ) )
			{
				broadcastMessage_warn( "System does not seem to contain devices other than device type \"" + std::to_string(device_type) + "\"!" );
			}
			else
			{
				broadcastMessage_fatal( "Got unknown error code. Aborting..." );
				failTest();
			}
		}
	}

	BOOST_AUTO_TEST_CASE(disableCpusByCommandLine)
	{
		const char * _params[] = {"foo", "--use_cpu=false"};
		disableSpecificDeviceTypeByCommandLine( _params, CL_DEVICE_TYPE_CPU);
	}

	BOOST_AUTO_TEST_CASE(disableGpusByCommandLine)
	{
		const char * _params[] = {"foo", "--use_gpu=false"};
		disableSpecificDeviceTypeByCommandLine( _params, CL_DEVICE_TYPE_GPU);
	}

BOOST_AUTO_TEST_SUITE_END()

void checkOnProperEnvironmentSettings()
{
	BOOST_REQUIRE_EQUAL(std::string("3"), std::string(getenv("GPU_DUMP_DEVICE_KERNEL")));
	BOOST_REQUIRE_NE(std::string(getenv("AMD_OCL_BUILD_OPTIONS_APPEND")).find("-save-temps"), -1);
}

BOOST_AUTO_TEST_CASE(dump_source_if_debugging)
{
	const char * _params[] = {"foo", "--log-level=debug"};
	meta::Inputparameters params(2, _params);
	hardware::System system(params);

	if(logger.beDebug()) {
		checkOnProperEnvironmentSettings();
	}
	else
	{
		broadcastMessage_fatal( "Something went wrong, logger not in debug mode..." );
		failTest();
	}
}
