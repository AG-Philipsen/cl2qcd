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

#include "../meta/util.hpp"
#include "../host_functionality/host_random.h"
#include "../hardware/code/spinors_staggered.hpp"
#include "../hardware/code/spinors.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE_FERMIONS_STAGGERED
#include <boost/test/unit_test.hpp>

//some functionality
#include "test_util.h"
#include "test_util_staggered.h"


void fill_sf_with_one(su3vec * sf_in, int size)
{
	for(int i = 0; i < size; ++i) {
		sf_in[i].e0 = hmc_complex_one;
		sf_in[i].e1 = hmc_complex_one;
		sf_in[i].e2 = hmc_complex_one;
	}
	return;
}

void fill_sf_with_random(su3vec * sf_in, int size, int seed)
{
	prng_init(seed);
	for(int i = 0; i < size; ++i) {
		sf_in[i].e0.re = prng_double();
		sf_in[i].e1.re = prng_double();
		sf_in[i].e2.re = prng_double();

		sf_in[i].e0.im = prng_double();
		sf_in[i].e1.im = prng_double();
		sf_in[i].e2.im = prng_double();
	}
	return;
}

void fill_sf_with_random(su3vec * sf_in, int size)
{
	fill_sf_with_random(sf_in, size, 123456);
}



void test_build(std::string inputfile)
{
	logger.info() << "build opencl_module_fermions_staggered";
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield_stagg cpu(&system);
	BOOST_MESSAGE("Test done");
}

void test_m_staggered(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "M_staggered";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield_stagg cpu(&system);
	
	//The following three lines are to be used to produce the ref_conf file needed to get the ref_value
	//---> Comment them out when the reference values have been obtained!
	/*
	print_gaugefield_to_textfile("ref_conf",&cpu,params);
	logger.info() << "Produced the ref_conf text file with the links for the ref. code. Returning...";
	return;
	// */

	cl_int err = CL_SUCCESS;
	const hardware::code::Fermions_staggered * device = cpu.get_device();
	su3vec * sf_in;
	su3vec * sf_out;

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_spinorfieldsize(params);

	sf_in = new su3vec[NUM_ELEMENTS_SF];
	sf_out = new su3vec[NUM_ELEMENTS_SF];

	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	else fill_sf_with_random(sf_in, NUM_ELEMENTS_SF);
	BOOST_REQUIRE(sf_in);
	
	//The following three lines are to be used to produce the ref_vec file needed to get the ref_value
	//---> Comment them out when the reference values have been obtained!
	/*
	print_staggeredfield_to_textfile("ref_vec",sf_in,params);
	logger.info() << "Produced the ref_vec text file with the staggered field for the ref. code. Returning...";
	return;
	 */

	const Plain<su3vec> in(NUM_ELEMENTS_SF, device->get_device());
	in.load(sf_in);
	const Plain<su3vec> out(NUM_ELEMENTS_SF, device->get_device());
	out.load(sf_in);
	const hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	auto spinor_code = device->get_device()->get_spinor_staggered_code();
		
	logger.info() << "|phi|^2:";
	hmc_float cpu_back;
	spinor_code->set_float_to_global_squarenorm_device(&in, &sqnorm);
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;
	logger.info() << "Run kernel";
	device->M_staggered_device(&in, &out,  cpu.get_gaugefield(), ARG_DEF);

	logger.info() << "result:";
	hmc_float cpu_res;
	spinor_code->set_float_to_global_squarenorm_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << std::setprecision(16) << cpu_res;
	
	logger.info() << "Clear buffers";
	delete[] sf_in;
	delete[] sf_out;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}


void test_DKS_eo(std::string inputfile)
{
	using namespace hardware::buffers;
	std::string kernelName = "D_KS_eo";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield_stagg cpu(&system);
	auto * device = cpu.get_device();
	
	//The following three lines are to be used to produce the ref_conf file needed to get the ref_value
	//---> Comment them out when the reference values have been obtained!
	/*
	print_gaugefield_to_textfile("ref_conf",&cpu,params);
	logger.info() << "Produced the ref_conf text file with the links for the ref. code. Returning...";
	//return;
	// */

	logger.info() << "Fill buffers...";
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	size_t NUM_ELEMENTS_SF_EO = hardware::code::get_eoprec_spinorfieldsize(params);
	su3vec * sf_in_eo;
	sf_in_eo = new su3vec[NUM_ELEMENTS_SF_EO];
	const SU3vec in_eo_even(NUM_ELEMENTS_SF_EO, device->get_device());
	const SU3vec out_eo(NUM_ELEMENTS_SF_EO, device->get_device());
	if(params.get_solver() == meta::Inputparameters::cg) fill_sf_with_one(sf_in_eo, NUM_ELEMENTS_SF_EO);
	else fill_sf_with_random(sf_in_eo, NUM_ELEMENTS_SF_EO);
	in_eo_even.load(sf_in_eo);
	
	//The following six lines are to be used to produce the ref_vec file needed to get the ref_value
	//---> Comment them out when the reference values have been obtained!
	/*
	if(params.get_read_multiple_configs())
	  print_staggeredfield_eo_to_textfile("ref_vec_odd",sf_in_eo,params);
	else
	  print_staggeredfield_eo_to_textfile("ref_vec_even",sf_in_eo,params);
	logger.info() << "Produced the ref_vec text file with the staggered field for the ref. code. Returning...";
	return;
	// */

	auto spinor_code = device->get_device()->get_spinor_staggered_code();

	logger.info() << "|phi|^2:";
	hmc_float cpu_back;
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in_eo_even, &sqnorm);
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;

	hmc_float cpu_res;
	if(params.get_read_multiple_configs()) {
		device->D_KS_eo_device( &in_eo_even, &out_eo, cpu.get_gaugefield(), EVEN);
	} else {
		device->D_KS_eo_device( &in_eo_even, &out_eo, cpu.get_gaugefield(), ODD);
	}
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&out_eo, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << "result:";
	logger.info() << cpu_res;

	logger.info() << "Clear buffers";
	delete[] sf_in_eo;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}



BOOST_AUTO_TEST_SUITE(BUILD)

BOOST_AUTO_TEST_CASE( BUILD_1 )
{
	test_build("/opencl_module_fermions_staggered_build_input_1");
}

BOOST_AUTO_TEST_CASE( BUILD_2 )
{
	test_build("/opencl_module_fermions_staggered_build_input_2");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( M_STAGGERED )

BOOST_AUTO_TEST_CASE( M_STAGGERED_1)
{
	test_m_staggered("/m_staggered_input_1");
}

BOOST_AUTO_TEST_CASE( M_STAGGERED_2)
{
	test_m_staggered("/m_staggered_input_2");
}

BOOST_AUTO_TEST_CASE( M_STAGGERED_3)
{
	test_m_staggered("/m_staggered_input_3");
}

BOOST_AUTO_TEST_CASE( M_STAGGERED_4)
{
	test_m_staggered("/m_staggered_input_4");
}

BOOST_AUTO_TEST_CASE( M_STAGGERED_5)
{
	test_m_staggered("/m_staggered_input_5");
}

BOOST_AUTO_TEST_CASE( M_STAGGERED_6)
{
	test_m_staggered("/m_staggered_input_6");
}

BOOST_AUTO_TEST_CASE( M_STAGGERED_7)
{
	test_m_staggered("/m_staggered_input_7");
}

BOOST_AUTO_TEST_CASE( M_STAGGERED_8)
{
	test_m_staggered("/m_staggered_input_8");
}

BOOST_AUTO_TEST_CASE( M_STAGGERED_9)
{
	test_m_staggered("/m_staggered_input_9");
}

BOOST_AUTO_TEST_CASE( M_STAGGERED_10)
{
	test_m_staggered("/m_staggered_input_10");
}

BOOST_AUTO_TEST_CASE( M_STAGGERED_11)
{
	test_m_staggered("/m_staggered_input_11");
}

BOOST_AUTO_TEST_CASE( M_STAGGERED_12)
{
	test_m_staggered("/m_staggered_input_12");
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( DKS_EO )

BOOST_AUTO_TEST_CASE( DKS_EO_1)
{
	test_DKS_eo("/dks_input_1");
}

BOOST_AUTO_TEST_CASE( DKS_EO_2)
{
	test_DKS_eo("/dks_input_2");
}

BOOST_AUTO_TEST_CASE( DKS_EO_3)
{
	test_DKS_eo("/dks_input_3");
}

BOOST_AUTO_TEST_CASE( DKS_EO_4)
{
	test_DKS_eo("/dks_input_4");
}

BOOST_AUTO_TEST_CASE( DKS_EO_5)
{
	test_DKS_eo("/dks_input_5");
}

BOOST_AUTO_TEST_CASE( DKS_EO_6)
{
	test_DKS_eo("/dks_input_6");
}

BOOST_AUTO_TEST_CASE( DKS_EO_7)
{
	test_DKS_eo("/dks_input_7");
}

BOOST_AUTO_TEST_CASE( DKS_EO_8)
{
	test_DKS_eo("/dks_input_8");
}

BOOST_AUTO_TEST_CASE( DKS_EO_9)
{
	test_DKS_eo("/dks_input_9");
}

BOOST_AUTO_TEST_CASE( DKS_EO_10)
{
	test_DKS_eo("/dks_input_10");
}

BOOST_AUTO_TEST_CASE( DKS_EO_11)
{
	test_DKS_eo("/dks_input_11");
}

BOOST_AUTO_TEST_CASE( DKS_EO_12)
{
	test_DKS_eo("/dks_input_12");
}

BOOST_AUTO_TEST_CASE( DKS_EO_13)
{
	test_DKS_eo("/dks_input_13");
}

BOOST_AUTO_TEST_CASE( DKS_EO_14)
{
	test_DKS_eo("/dks_input_14");
}

BOOST_AUTO_TEST_CASE( DKS_EO_15)
{
	test_DKS_eo("/dks_input_15");
}

BOOST_AUTO_TEST_CASE( DKS_EO_16)
{
	test_DKS_eo("/dks_input_16");
}

BOOST_AUTO_TEST_SUITE_END()

