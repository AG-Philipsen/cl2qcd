/**
 * Copyright (c) 2016 Christopher Pinke
 * Copyright (c) 2018,2021 Alessandro Sciarra
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.
 */

#include "hardwareTestUtilities.hpp"

#include <boost/program_options.hpp>
#include <boost/test/unit_test.hpp>

static std::tuple<bool, bool, bool> parseBoostRuntimeArguments()
{
    static bool parsed   = false;
    static bool useGpu   = false;
    static bool useCpu   = false;
    static bool useRec12 = false;
    auto argc            = boost::unit_test::framework::master_test_suite().argc;
    auto argv            = boost::unit_test::framework::master_test_suite().argv;
    if (parsed == false && argc > 1) {  // argv[0] is the executable name
        namespace po = boost::program_options;
        logger.info() << "Found " << argc << " runtime user command line arguments";
        po::options_description commandLineOptionsDescription("Test command line options");
        // clang-format off
        commandLineOptionsDescription.add_options()
        ("useGPU", po::value<bool>(&useGpu)->default_value(useGpu), "Whether to use GPUs in tests. If set to true and at least one GPU is found, GPUs only are used.")
        ("useCPU", po::value<bool>(&useCpu)->default_value(useCpu), "Whether to use CPUs in tests.")
        ("useReconstruct12", po::value<bool>(&useRec12)->default_value(useRec12), "Whether to use the reconstruction technique for the gaugefield in tests.");
        // clang-format on
        po::variables_map vm;
        try {
            po::store(po::parse_command_line(argc, argv, commandLineOptionsDescription), vm);
            boost::program_options::notify(vm);
        } catch (po::error& e) {
            logger.error() << e.what() << "\n\n" << commandLineOptionsDescription << "\n";
            exit(1);
        }
        parsed = true;
    }
    return std::make_tuple(useCpu, useGpu, useRec12);
}

bool checkBoostRuntimeArgumentsForGpuUsage()
{
    return std::get<1>(parseBoostRuntimeArguments());
}
bool checkBoostRuntimeArgumentsForRec12Usage()
{
    return std::get<2>(parseBoostRuntimeArguments());
}

void broadcastMessage_warn(const std::string message)
{
    logger.warn() << message;
    BOOST_TEST_MESSAGE(message);
}

void broadcastMessage_error(const std::string message)
{
    logger.error() << message;
    BOOST_TEST_MESSAGE(message);
}

void broadcastMessage_fatal(const std::string message)
{
    logger.fatal() << message;
    BOOST_TEST_MESSAGE(message);
}

void failTest()
{
    BOOST_ERROR("Test failure!");
}

void atLeastOneDeviceMustExistForSanityOfSystem(const hardware::System* system)
{
    BOOST_REQUIRE_GE(system->get_devices().size(), 1);
}

bool checkIfNoOpenCLDevicesWereFound(const hardware::OpenclException exception)
{
    return exception.errorCode == -1;
}

void endTestAsNoDevicesWereFound()
{
    broadcastMessage_warn("System does not seem to contain devices of desired type!");
    broadcastMessage_warn("Exiting...");
    exit(0);
}

void endTestBecauseOfUnknownError()
{
    broadcastMessage_fatal("Got unknown error code. Aborting...");
    failTest();
}

void handleExceptionInTest(hardware::OpenclException& exception)
{
    if (checkIfNoOpenCLDevicesWereFound(exception)) {
        endTestAsNoDevicesWereFound();
    } else {
        endTestBecauseOfUnknownError();
    }
}
