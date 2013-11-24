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

#include <iostream>
// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE staple_test
#include <boost/test/unit_test.hpp>

void printWorkaroundInfo(){
	std::cout << "If this test fails, a workaround can be used that performs the test with std args." << std::endl;
	std::cout << "To enable this, set" << std::endl;
	std::cout << "\tBOOST_FAILURE_WORKAROUND" << std::endl;
	std::cout << "to 1 in \"test_util.h\"" << std::endl;
}

BOOST_AUTO_TEST_CASE( BOOST_ARGC_1 )
{
	int param_expect = 1;
	std::cout << "Testing getting args from boost master_test_suite..." << std::endl;
	std::cout << "Getting number of arguments, this should be equal to " << param_expect << "!" << std::endl;
	printWorkaroundInfo();
	int num_par = 	boost::unit_test::framework::master_test_suite().argc;
	BOOST_REQUIRE_EQUAL(num_par, param_expect);
	BOOST_MESSAGE("Test done");
}

BOOST_AUTO_TEST_CASE( BOOST_ARGC_2 )
{
	int param_expect = 2;
	std::cout << "Testing getting args from boost master_test_suite..." << std::endl;
	std::cout << "Getting number of arguments, this should be equal to " << param_expect << "!" << std::endl;
	printWorkaroundInfo();
	int num_par = 	boost::unit_test::framework::master_test_suite().argc;
	BOOST_REQUIRE_EQUAL(num_par, param_expect);
	BOOST_MESSAGE("Test done");
}

BOOST_AUTO_TEST_CASE( BOOST_ARGV )
{
	int param_expect = 4;
	int num_par = 	boost::unit_test::framework::master_test_suite().argc;
	printWorkaroundInfo();
	BOOST_REQUIRE_EQUAL(param_expect, num_par);
	if(num_par < param_expect){
		std::cout << "Got only " << num_par << " Inputparameters, expected " << param_expect << "! Use inputfile values instead!"<< std::endl;
	} else if(num_par > param_expect){
		std::cout << "Got " << num_par << " Inputparameters, expected only " << param_expect << "! Use only the first " << param_expect << " values" << std::endl;
	}
	
	//get input file location that has been passed as an argument
	std::string inputfile_location = boost::unit_test::framework::master_test_suite().argv[1];
	std::cout << "inputfile used: " << inputfile_location<< std::endl;
	//get use_gpu = true/false that has been passed as an argument
	std::string gpu_opt =  boost::unit_test::framework::master_test_suite().argv[2];
	std::cout << "GPU usage: " << gpu_opt << std::endl;
	//get use_rec12 = true/false that has been passed as an argument
	std::string rec12_opt =  boost::unit_test::framework::master_test_suite().argv[3];
	std::cout << "rec12 usage: " << rec12_opt << std::endl;
}
