/*
 * Copyright (c) 2012,2014 Christopher Pinke
 * Copyright (c) 2013 Matthias Bach
 * Copyright (c) 2018 Alessandro Sciarra
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

#include "../host_functionality/logger.hpp"

#include <iostream>

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE staple_test
#include <boost/test/unit_test.hpp>
#include <boost/version.hpp>

/*
 * TODO: Since version 1.60 of boost,
 *          https://www.boost.org/doc/libs/1_60_0/libs/test/doc/html/boost_test/change_log.html
 *       the boost command line parameters have been forced to be separated from the user program ones using a '--'
 *       on the command line. Since version 1.60 Such a '--' is not figuring in the argc/argv of the master_test_suite
 *       of boost while in previous version it does. Therefore, there might be an inconsistency between different
 *       installations. This has been taken into account here using the version.hpp header. Note that at the moment
 *       we require version 1.59 of boost in order to use smart filtering in the hardware::System test, but everywhere
 *       we already prepared the code base to use the new '--' convention introduced with boost 1.60. The next step will
 *       be to require such a version, cleaning up the code here.
 */

void printInfo(int expectedNumberOfParameters)
{
    logger.info() << "Testing passing of arguments from boost master_test_suite...";
    logger.info() << "Getting number of arguments, this should be equal to " << expectedNumberOfParameters << "!";
    logger.info() << "If this test fails, it might indicate a problem of boost with the used compiler...";
}

void checkArgc(int expectedNumberOfParameters)
{
    printInfo(expectedNumberOfParameters);
    int numberOfParametersFromBoost = boost::unit_test::framework::master_test_suite().argc;
    BOOST_REQUIRE_EQUAL(numberOfParametersFromBoost, expectedNumberOfParameters);
}

BOOST_AUTO_TEST_SUITE(BOOST_ARGUMENTS)

    BOOST_AUTO_TEST_CASE(BOOST_ARGC_1)
    {
        // Run as: ./boost_test --run_test=BOOST_ARGUMENTS/BOOST_ARGC_1
        int expectedNumberOfParameters = 1;
        checkArgc(expectedNumberOfParameters);
    }

    BOOST_AUTO_TEST_CASE(BOOST_ARGC_2)
    {
        // Run as: ./boost_test --run_test=BOOST_ARGUMENTS/BOOST_ARGC_2 -- "firstParam"
        int expectedNumberOfParameters = 2;
        if (BOOST_VERSION < 106000)
            expectedNumberOfParameters = 3;  // boost args still count (plus '--')!
        checkArgc(expectedNumberOfParameters);
    }

    void checkArgv(int position, std::string expectedContent)
    {
        std::string argument = boost::unit_test::framework::master_test_suite().argv[position];
        BOOST_REQUIRE_EQUAL(argument, expectedContent);
    }

    BOOST_AUTO_TEST_CASE(BOOST_ARGV)
    {
        // Run as:
        // ./boost_test --run_test=BOOST_ARGUMENTS/BOOST_ARGC_2 -- "firstArgument" "secondArgument" "thirdArgument"
        int expectedNumberOfParameters = 4, argIndex = 1;

        if (BOOST_VERSION < 106000) {
            expectedNumberOfParameters = 5;  // boost args still count (plus '--')!
            checkArgv(argIndex++, "--");
        }
        checkArgc(expectedNumberOfParameters);

        checkArgv(argIndex++, "firstArgument");
        checkArgv(argIndex++, "secondArgument");
        checkArgv(argIndex++, "thirdArgument");
    }

BOOST_AUTO_TEST_SUITE_END()
