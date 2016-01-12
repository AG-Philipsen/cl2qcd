/**
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

#include "hardwareTestUtilities.hpp"
#include <boost/test/unit_test.hpp>

std::pair<bool, bool> checkForBoostRuntimeArguments()
{
	bool useGpu = false;
	bool useRec12 = false;
	int num_par = boost::unit_test::framework::master_test_suite().argc;
	if(num_par != 0)
	{
		logger.info() << "Found " << num_par << " runtime arguments, checking for gpu and rec12 options...";
		for (int i = 1; i< num_par; i++)
		{
			std::string currentArgument = boost::unit_test::framework::master_test_suite().argv[i];
			if (currentArgument.find("--use_gpu") != std::string::npos)
			{
				if(currentArgument.find("true") != std::string::npos)
				{
					useGpu = true;
				}
			}
			if (currentArgument.find("--use_rec12") != std::string::npos)
			{
				if(currentArgument.find("true") != std::string::npos)
				{
					useRec12 = true;
				}
			}
		}
	}
	return std::pair<bool, bool>{useGpu, useRec12};
}

bool checkBoostRuntimeArgumentsForGpuUsage()
{
	return checkForBoostRuntimeArguments().first;
}
bool checkBoostRuntimeArgumentsForRec12Usage()
{
	return checkForBoostRuntimeArguments().second;
}

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
	BOOST_CHECK_EQUAL(true, false);
}

void atLeastOneDeviceMustExistForSanityOfSystem(const hardware::System * system)
{
	BOOST_REQUIRE_GE(system->get_devices().size(), 1);
}

bool checkIfNoOpenCLDevicesWereFound( const hardware::OpenclException exception)
{
	return exception.errorCode == -1;
}

void endTestAsNoDevicesWereFound()
{
	broadcastMessage_warn( "System does not seem to contain devices of desired type!" );
	broadcastMessage_warn( "Exiting..." );
	exit(0);
}

void endTestBecauseOfUnknownError()
{
	broadcastMessage_fatal( "Got unknown error code. Aborting..." );
	failTest();
}

void handleExceptionInTest(hardware::OpenclException & exception)
{
	if ( checkIfNoOpenCLDevicesWereFound( exception ) )
	{
		endTestAsNoDevicesWereFound();
	}
	else
	{
		endTestBecauseOfUnknownError();
	}
}

