#include <iostream>
// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE staple_test
#include <boost/test/unit_test.hpp>

BOOST_AUTO_TEST_CASE( BOOST_ARGC_1 )
{
	int param_expect = 1;
	std::cout << "Testing getting args from boost master_test_suite..." << std::endl;
	std::cout << "Getting number of arguments, this should be equal to " << param_expect << "!" << std::endl;
	int num_par = 	boost::unit_test::framework::master_test_suite().argc;
	BOOST_REQUIRE_EQUAL(num_par, param_expect);
	BOOST_MESSAGE("Test done");
}

BOOST_AUTO_TEST_CASE( BOOST_ARGC_2 )
{
	int param_expect = 2;
	std::cout << "Testing getting args from boost master_test_suite..." << std::endl;
	std::cout << "Getting number of arguments, this should be equal to " << param_expect << "!" << std::endl;
	int num_par = 	boost::unit_test::framework::master_test_suite().argc;
	BOOST_REQUIRE_EQUAL(num_par, param_expect);
	BOOST_MESSAGE("Test done");
}

BOOST_AUTO_TEST_CASE( BOOST_ARGV )
{
	int param_expect = 4;
	int num_par = 	boost::unit_test::framework::master_test_suite().argc;
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
