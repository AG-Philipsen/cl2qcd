/*
 * Copyright 2012, 2013, 2014 Lars Zeidlewicz, Christopher Pinke,
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

#include <iostream>
// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE staple_test
#include <boost/test/unit_test.hpp>

void printInfo(int expectedNumberOfParameters){
	std::cout << "Testing passing of arguments from boost master_test_suite..." << std::endl;
	std::cout << "Getting number of arguments, this should be equal to " << expectedNumberOfParameters << "!" << std::endl;
	std::cout << "If this test fails, it might indicate a problem of boost with the used compiler..." << std::endl;
}

void checkArgc(int expectedNumberOfParameters)
{
		printInfo(expectedNumberOfParameters);
		int numberOfParametersFromBoost = 	boost::unit_test::framework::master_test_suite().argc;
		BOOST_REQUIRE_EQUAL(numberOfParametersFromBoost, expectedNumberOfParameters);
}

BOOST_AUTO_TEST_SUITE ( BOOST_ARGUMENTS )

	BOOST_AUTO_TEST_CASE( BOOST_ARGC_1 )
	{
		int expectedNumberOfParameters = 1;
		checkArgc(expectedNumberOfParameters);
	}
	
	BOOST_AUTO_TEST_CASE( BOOST_ARGC_2 )
	{
		int expectedNumberOfParameters = 2;
		checkArgc(expectedNumberOfParameters);
	}

	void checkArgv(int position, std::string expectedContent)
	{
		std::string argument = boost::unit_test::framework::master_test_suite().argv[position];
		BOOST_REQUIRE_EQUAL(argument, expectedContent);
	}
	
	BOOST_AUTO_TEST_CASE( BOOST_ARGV )
	{
		int expectedNumberOfParameters = 4;
		checkArgc(expectedNumberOfParameters);
		
		checkArgv(1, "firstArgument");
		checkArgv(2, "secondArgument");
		checkArgv(3, "thirdArgument");
	}
	
BOOST_AUTO_TEST_SUITE_END()