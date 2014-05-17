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

#include "testUtilities.hpp"

#include "../meta/util.hpp"
#include "../host_functionality/host_random.h"
#include "../physics/prng.hpp"
#include "../hardware/device.hpp"
#include "../hardware/code/spinors.hpp"
#include "../hardware/code/complex.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE_COMPLEX
#include <boost/test/unit_test.hpp>

//some functionality
#include "../../tests/test_util.h"

void test_build(std::string inputfile)
{
	logger.info() << "build opencl_module_complex";
	logger.info() << "Init device";
	meta::Inputparameters params = createParameters("complex/" + inputfile);
	hardware::System system(params);
	for(auto device: system.get_devices()) {
		device->get_complex_code();
	}
	BOOST_MESSAGE("Test done");
}

void test_cplx(std::string inputfile, int switcher, bool hard=false)
{
  //switcher chooses between product, ratio, sum, subtraction and convert
	using namespace hardware::buffers;

	std::string kernelName;
	if (switcher == 0)
	  kernelName = "product";
	else if(switcher == 1)
	  kernelName = "ratio";
	else if(switcher == 2)
	  kernelName = "sum";
	else if(switcher == 3)
	  kernelName = "subtraction";
	else if (switcher == 4)
	  kernelName = "convert";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = createParameters("complex/" + inputfile);
	hardware::System system(params);
	auto * device = system.get_devices().at(0)->get_complex_code();

	logger.info() << "Fill buffers...";
	hardware::buffers::Plain<hmc_complex> sqnorm(1, device->get_device());
	hardware::buffers::Plain<hmc_complex> alpha(1, device->get_device());
	hardware::buffers::Plain<hmc_complex> beta(1, device->get_device());

	hmc_complex alpha_host = {params.get_beta(), params.get_rho()};
	logger.info() << "Use alpha = (" << alpha_host.re << ","<< alpha_host.im <<")";
	hmc_complex beta_host = {params.get_kappa(), params.get_mu()};
	logger.info() << "Use beta = (" << beta_host.re << ","<< beta_host.im <<")";

	alpha.load(&alpha_host);
	beta.load(&beta_host);

	logger.info() << "Run kernel";
	if(switcher == 0){
	  device->set_complex_to_product_device(&alpha, &beta, &sqnorm);
	  if(hard)
	    device->set_complex_to_product_device(&alpha, &sqnorm, &sqnorm);
	}
	else if (switcher ==1){
	  device->set_complex_to_ratio_device(&alpha, &beta, &sqnorm);
	  if(hard)
	    device->set_complex_to_ratio_device(&sqnorm, &beta, &sqnorm);
	}
	else if (switcher ==2){
	  device->set_complex_to_sum_device(&alpha, &beta, &sqnorm);
	  if(hard)
	    device->set_complex_to_sum_device(&alpha, &sqnorm, &sqnorm);
	}
	else if (switcher ==3){
	  device->set_complex_to_difference_device(&alpha, &beta, &sqnorm);
	  if(hard)
	    device->set_complex_to_difference_device(&sqnorm, &beta, &sqnorm);
	}
	if(switcher == 4){
	  hardware::buffers::Plain<hmc_float> gamma(1, device->get_device());
	  hmc_float tmp = (params.get_beta());
	  gamma.load(&tmp);
	  device->set_complex_to_float_device(&gamma, &sqnorm);
	}
	logger.info() << "result:";
	hmc_float cpu_res;
	hmc_complex tmp;
	sqnorm.dump(&tmp);
	cpu_res = tmp.re + tmp.im;
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}


BOOST_AUTO_TEST_SUITE(BUILD)

BOOST_AUTO_TEST_CASE( BUILD_1 )
{
	test_build("complex_build_input_1");
}

BOOST_AUTO_TEST_CASE( BUILD_2 )
{
	test_build("complex_build_input_2");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CPLX_PRODUCT)

BOOST_AUTO_TEST_CASE( CPLX_PRODUCT_1 )
{
  test_cplx("/cplx_product_input_1", 0);
}

BOOST_AUTO_TEST_CASE( CPLX_PRODUCT_2 )
{
  test_cplx("/cplx_product_input_2", 0);
}

BOOST_AUTO_TEST_CASE( CPLX_PRODUCT_3 )
{
  test_cplx("/cplx_product_input_3", 0);
}

BOOST_AUTO_TEST_CASE( CPLX_PRODUCT_4 )
{
  test_cplx("/cplx_product_input_4", 0);
}

BOOST_AUTO_TEST_CASE( CPLX_PRODUCT_5 )
{
  test_cplx("/cplx_product_input_5", 0);
}

BOOST_AUTO_TEST_CASE( CPLX_PRODUCT_6 )
{
  test_cplx("/cplx_product_input_6", 0);
}

BOOST_AUTO_TEST_CASE( CPLX_PRODUCT_7 )
{
  test_cplx("/cplx_product_input_7", 0, true);
}

BOOST_AUTO_TEST_CASE( CPLX_PRODUCT_8 )
{
  test_cplx("/cplx_product_input_8", 0, true);
}

BOOST_AUTO_TEST_CASE( CPLX_PRODUCT_9 )
{
  test_cplx("/cplx_product_input_9", 0, true);
}

BOOST_AUTO_TEST_CASE( CPLX_PRODUCT_10 )
{
  test_cplx("/cplx_product_input_10", 0, true);
}

BOOST_AUTO_TEST_CASE( CPLX_PRODUCT_11 )
{
  test_cplx("/cplx_product_input_11", 0, true);
}

BOOST_AUTO_TEST_CASE( CPLX_PRODUCT_12 )
{
  test_cplx("/cplx_product_input_12", 0, true);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CPLX_RATIO)

BOOST_AUTO_TEST_CASE( CPLX_RATIO_1 )
{
  test_cplx("/cplx_ratio_input_1", 1);
}

BOOST_AUTO_TEST_CASE( CPLX_RATIO_2 )
{
  test_cplx("/cplx_ratio_input_2", 1);
}

BOOST_AUTO_TEST_CASE( CPLX_RATIO_3 )
{
  test_cplx("/cplx_ratio_input_3", 1);
}

BOOST_AUTO_TEST_CASE( CPLX_RATIO_4 )
{
  test_cplx("/cplx_ratio_input_4", 1);
}

BOOST_AUTO_TEST_CASE( CPLX_RATIO_5 )
{
  test_cplx("/cplx_ratio_input_5", 1, true);
}

BOOST_AUTO_TEST_CASE( CPLX_RATIO_6 )
{
  test_cplx("/cplx_ratio_input_6", 1, true);
}

BOOST_AUTO_TEST_CASE( CPLX_RATIO_7 )
{
  test_cplx("/cplx_ratio_input_7", 1, true);
}

BOOST_AUTO_TEST_CASE( CPLX_RATIO_8 )
{
  test_cplx("/cplx_ratio_input_8", 1, true);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CPLX_SUM)

BOOST_AUTO_TEST_CASE( CPLX_SUM_1 )
{
  test_cplx("/cplx_sum_input_1", 2);
}

BOOST_AUTO_TEST_CASE( CPLX_SUM_2 )
{
  test_cplx("/cplx_sum_input_2", 2);
}

BOOST_AUTO_TEST_CASE( CPLX_SUM_3 )
{
  test_cplx("/cplx_sum_input_3", 2);
}

BOOST_AUTO_TEST_CASE( CPLX_SUM_4 )
{
  test_cplx("/cplx_sum_input_4", 2, true);
}

BOOST_AUTO_TEST_CASE( CPLX_SUM_5 )
{
  test_cplx("/cplx_sum_input_5", 2, true);
}

BOOST_AUTO_TEST_CASE( CPLX_SUM_6 )
{
  test_cplx("/cplx_sum_input_6", 2, true);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CPLX_DIFFERENCE)

BOOST_AUTO_TEST_CASE( CPLX_DIFFERENCE_1 )
{
  test_cplx("/cplx_difference_input_1", 3);
}

BOOST_AUTO_TEST_CASE( CPLX_DIFFERENCE_2 )
{
  test_cplx("/cplx_difference_input_2", 3);
}

BOOST_AUTO_TEST_CASE( CPLX_DIFFERENCE_3 )
{
  test_cplx("/cplx_difference_input_3", 3);
}

BOOST_AUTO_TEST_CASE( CPLX_DIFFERENCE_4 )
{
  test_cplx("/cplx_difference_input_4", 3, true);
}

BOOST_AUTO_TEST_CASE( CPLX_DIFFERENCE_5 )
{
  test_cplx("/cplx_difference_input_5", 3, true);
}

BOOST_AUTO_TEST_CASE( CPLX_DIFFERENCE_6 )
{
  test_cplx("/cplx_difference_input_6", 3, true);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CPLX_CONVERT)

BOOST_AUTO_TEST_CASE( CPLX_CONVERT_1 )
{
  test_cplx("/cplx_convert_input_1", 4);
}

BOOST_AUTO_TEST_CASE( CPLX_CONVERT_2 )
{
  test_cplx("/cplx_convert_input_2", 4);
}

BOOST_AUTO_TEST_SUITE_END()
