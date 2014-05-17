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
#include "../physics/prng.hpp"
#include "../hardware/device.hpp"
#include "../hardware/code/real.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE_REAL
#include <boost/test/unit_test.hpp>

//some functionality
#include "../../tests/test_util.h"

void test_build(std::string inputfile)
{
	logger.info() << "build opencl_module_real";
	logger.info() << "Init device";
	meta::Inputparameters params = createParameters("real/" + inputfile);
	hardware::System system(params);
	for(auto device: system.get_devices()) {
		device->get_real_code();
	}
	BOOST_MESSAGE("Test done");
}

void test_access_element(std::string inputfile, bool get=true)
{
	using namespace hardware::buffers;
	
	std::string kernelName;
	if (get)
	  kernelName = "get_vector_element";
	else
	  kernelName = "set_vector_element";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = createParameters("real/" + inputfile);
	hardware::System system(params);
	auto * device = system.get_devices().at(0)->get_real_code();
	
	logger.info() << "Fill buffers...";
	hardware::buffers::Plain<hmc_float> scalar_buf(1, device->get_device());
	hardware::buffers::Plain<hmc_float> vector_buf(4, device->get_device());
	std::vector<hmc_float> vector_host;
	hmc_float scalar_host;
	
	vector_host.push_back(params.get_beta());
	vector_host.push_back(params.get_kappa());
	vector_host.push_back(params.get_rho());
	vector_host.push_back(params.get_mu());
	vector_buf.load(&vector_host[0]);
	if(get){
	  logger.info() << "Using:";
	  for(uint i=0; i<vector_host.size(); i++)
	    logger.info() << "vector_buf[" << i << "] = " << vector_host[i];
	}else{
	  scalar_host = params.get_beta() + params.get_kappa() + params.get_rho() + params.get_mu();
	  logger.info() << "Using scalar = " << scalar_host;
	  scalar_buf.load(&scalar_host);
	}
	
	logger.info() << "Run kernel";
	
	if(get){
	  for(uint i=0; i<vector_host.size(); i++){
	    device->set_real_to_vector_element_device(&vector_buf, i, &scalar_buf);
	    hmc_float cpu_res;
	    scalar_buf.dump(&cpu_res);
	    logger.info() << "element " << i << " = " << cpu_res;
	    BOOST_REQUIRE_CLOSE(cpu_res, vector_host[i], 1.e-8);
	  }
	}else{
	  for(uint i=0; i<vector_host.size(); i++){
	    device->set_vector_element_to_real_device(&scalar_buf, i, &vector_buf);
	    std::vector<hmc_float> cpu_res(4);
	    vector_buf.dump(&cpu_res[0]);
	    logger.info() << "element " << i << " = " << cpu_res[i];
	    BOOST_REQUIRE_CLOSE(cpu_res[i], scalar_host, 1.e-8);
	  }
	}
	
	BOOST_MESSAGE("Test done");
}

void test_base_operations(std::string inputfile, int switcher, bool hard=false)
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
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = createParameters("real/" + inputfile);
	hardware::System system(params);
	auto * device = system.get_devices().at(0)->get_real_code();

	logger.info() << "Fill buffers...";
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	hardware::buffers::Plain<hmc_float> alpha(1, device->get_device());
	hardware::buffers::Plain<hmc_float> beta(1, device->get_device());

	hmc_float alpha_host = params.get_beta();
	logger.info() << "Use alpha = " << alpha_host;
	hmc_float beta_host = params.get_kappa();
	logger.info() << "Use beta = " << beta_host;

	alpha.load(&alpha_host);
	beta.load(&beta_host);

	logger.info() << "Run kernel";
	if(switcher == 0){
	  device->set_real_to_product_device(&alpha, &beta, &sqnorm);
	  if(hard)
	    device->set_real_to_product_device(&alpha, &sqnorm, &sqnorm);
	}
	else if (switcher ==1){
	  device->set_real_to_ratio_device(&alpha, &beta, &sqnorm);
	  if(hard)
	    device->set_real_to_ratio_device(&sqnorm, &beta, &sqnorm);
	}
	else if (switcher ==2){
	  device->set_real_to_sum_device(&alpha, &beta, &sqnorm);
	  if(hard)
	    device->set_real_to_sum_device(&alpha, &sqnorm, &sqnorm);
	}
	else if (switcher ==3){
	  device->set_real_to_difference_device(&alpha, &beta, &sqnorm);
	  if(hard)
	    device->set_real_to_difference_device(&sqnorm, &beta, &sqnorm);
	}
	logger.info() << "result:";
	hmc_float cpu_res;
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}


void test_update(std::string inputfile, int switcher)
{
  //switcher chooses between alpha, beta and zeta update
	using namespace hardware::buffers;

	std::string kernelName;
	if (switcher == 0)
	  kernelName = "update_alpha_cgm";
	else if(switcher == 1)
	  kernelName = "update_beta_cgm";
	else if(switcher == 2)
	  kernelName = "update_zeta_cgm";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = createParameters("real/" + inputfile);
	hardware::System system(params);
	auto * device = system.get_devices().at(0)->get_real_code();

	logger.info() << "Fill buffers...";
	hardware::buffers::Plain<hmc_float> out(1, device->get_device());
	hardware::buffers::Plain<hmc_float> variable_1(1, device->get_device());
	hardware::buffers::Plain<hmc_float> variable_2(1, device->get_device());
	hardware::buffers::Plain<hmc_float> variable_3(1, device->get_device());
	hardware::buffers::Plain<hmc_float> variable_4(1, device->get_device());
	hardware::buffers::Plain<hmc_float> variable_5(1, device->get_device());
	hardware::buffers::Plain<hmc_float> variable_6(1, device->get_device());

	hmc_float variable_1_host = params.get_beta();
	hmc_float variable_2_host = params.get_rho();
	hmc_float variable_3_host = params.get_kappa();
	hmc_float variable_4_host = params.get_mu();
	hmc_float variable_5_host = params.get_mass(); 
	hmc_float variable_6_host = params.get_approx_lower();  
	
	logger.info() << "Use variable_1 = " << variable_1_host;
	logger.info() << "Use variable_2 = " << variable_2_host;
	logger.info() << "Use variable_3 = " << variable_3_host;
	if(switcher == 0 || switcher == 2){
	  logger.info() << "Use variable_4 = " << variable_4_host;
	  logger.info() << "Use variable_5 = " << variable_5_host;
	}
	if(switcher == 2)
	  logger.info() << "Use variable_6 = " << variable_6_host;

	variable_1.load(&variable_1_host);
	variable_2.load(&variable_2_host);
	variable_3.load(&variable_3_host);
	variable_4.load(&variable_4_host);
	variable_5.load(&variable_5_host);
	variable_6.load(&variable_6_host);

	logger.info() << "Run kernel";
	if(switcher == 0)
	  device->update_alpha_cgm_device(&variable_1, &variable_2, &variable_3, &variable_4, &variable_5, 1, &out);
	else if (switcher ==1)
	  device->update_beta_cgm_device(&variable_1, &variable_2, &variable_3, 1, &out);
	else if (switcher ==2)
	  device->update_zeta_cgm_device(&variable_1, &variable_2, &variable_3, &variable_4, &variable_5, &variable_6, 1, &out);
	logger.info() << "result:";
	hmc_float cpu_res;
	out.dump(&cpu_res);
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}


BOOST_AUTO_TEST_SUITE(BUILD)

BOOST_AUTO_TEST_CASE( BUILD_1 )
{
	test_build("/real_build_input_1"); 
}

BOOST_AUTO_TEST_CASE( BUILD_2 )
{
	test_build("/real_build_input_2");
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(ACCESS_ELEMENT)

BOOST_AUTO_TEST_CASE( ACCESS_ELEMENT_1 )
{
	test_access_element("/real_access_element_vector_input_1", true);
}

BOOST_AUTO_TEST_CASE( ACCESS_ELEMENT_2 )
{
	test_access_element("/real_access_element_vector_input_1", false);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(REAL_PRODUCT)

BOOST_AUTO_TEST_CASE( REAL_PRODUCT_1 )
{
  test_base_operations("/real_product_input_1", 0);
}

BOOST_AUTO_TEST_CASE( REAL_PRODUCT_2 )
{
  test_base_operations("/real_product_input_2", 0);
}

BOOST_AUTO_TEST_CASE( REAL_PRODUCT_3 )
{
  test_base_operations("/real_product_input_3", 0);
}

BOOST_AUTO_TEST_CASE( REAL_PRODUCT_4 )
{
  test_base_operations("/real_product_input_4", 0);
}

BOOST_AUTO_TEST_CASE( REAL_PRODUCT_5 )
{
  test_base_operations("/real_product_input_5", 0);
}

BOOST_AUTO_TEST_CASE( REAL_PRODUCT_6 )
{
  test_base_operations("/real_product_input_6", 0);
}

BOOST_AUTO_TEST_CASE( REAL_PRODUCT_7 )
{
  test_base_operations("/real_product_input_7", 0, true);
}

BOOST_AUTO_TEST_CASE( REAL_PRODUCT_8 )
{
  test_base_operations("/real_product_input_8", 0, true);
}

BOOST_AUTO_TEST_CASE( REAL_PRODUCT_9 )
{
  test_base_operations("/real_product_input_9", 0, true);
}

BOOST_AUTO_TEST_CASE( REAL_PRODUCT_10 )
{
  test_base_operations("/real_product_input_10", 0, true);
}

BOOST_AUTO_TEST_CASE( REAL_PRODUCT_11 )
{
  test_base_operations("/real_product_input_11", 0, true);
}

BOOST_AUTO_TEST_CASE( REAL_PRODUCT_12 )
{
  test_base_operations("/real_product_input_12", 0, true);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(REAL_RATIO)

BOOST_AUTO_TEST_CASE( REAL_RATIO_1 )
{
  test_base_operations("/real_ratio_input_1", 1);
}

BOOST_AUTO_TEST_CASE( REAL_RATIO_2 )
{
  test_base_operations("/real_ratio_input_2", 1);
}

BOOST_AUTO_TEST_CASE( REAL_RATIO_3 )
{
  test_base_operations("/real_ratio_input_3", 1);
}

BOOST_AUTO_TEST_CASE( REAL_RATIO_4 )
{
  test_base_operations("/real_ratio_input_4", 1);
}

BOOST_AUTO_TEST_CASE( REAL_RATIO_5 )
{
  test_base_operations("/real_ratio_input_5", 1, true);
}

BOOST_AUTO_TEST_CASE( REAL_RATIO_6 )
{
  test_base_operations("/real_ratio_input_6", 1, true);
}

BOOST_AUTO_TEST_CASE( REAL_RATIO_7 )
{
  test_base_operations("/real_ratio_input_7", 1, true);
}

BOOST_AUTO_TEST_CASE( REAL_RATIO_8 )
{
  test_base_operations("/real_ratio_input_8", 1, true);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(REAL_SUM)

BOOST_AUTO_TEST_CASE( REAL_SUM_1 )
{
  test_base_operations("/real_sum_input_1", 2);
}

BOOST_AUTO_TEST_CASE( REAL_SUM_2 )
{
  test_base_operations("/real_sum_input_2", 2);
}

BOOST_AUTO_TEST_CASE( REAL_SUM_3 )
{
  test_base_operations("/real_sum_input_3", 2);
}

BOOST_AUTO_TEST_CASE( REAL_SUM_4 )
{
  test_base_operations("/real_sum_input_4", 2, true);
}

BOOST_AUTO_TEST_CASE( REAL_SUM_5 )
{
  test_base_operations("/real_sum_input_5", 2, true);
}

BOOST_AUTO_TEST_CASE( REAL_SUM_6 )
{
  test_base_operations("/real_sum_input_6", 2, true);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(REAL_DIFFERENCE)

BOOST_AUTO_TEST_CASE( REAL_DIFFERENCE_1 )
{
  test_base_operations("/real_difference_input_1", 3);
}

BOOST_AUTO_TEST_CASE( REAL_DIFFERENCE_2 )
{
  test_base_operations("/real_difference_input_2", 3);
}

BOOST_AUTO_TEST_CASE( REAL_DIFFERENCE_3 )
{
  test_base_operations("/real_difference_input_3", 3);
}

BOOST_AUTO_TEST_CASE( REAL_DIFFERENCE_4 )
{
  test_base_operations("/real_difference_input_4", 3, true);
}

BOOST_AUTO_TEST_CASE( REAL_DIFFERENCE_5 )
{
  test_base_operations("/real_difference_input_5", 3, true);
}

BOOST_AUTO_TEST_CASE( REAL_DIFFERENCE_6 )
{
  test_base_operations("/real_difference_input_6", 3, true);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(REAL_UPDATE)

BOOST_AUTO_TEST_CASE( ALPHA_1 )
{
  test_update("/real_update_alpha_input_1", 0);
}

BOOST_AUTO_TEST_CASE( BET_1 )
{
  test_update("/real_update_beta_input_1", 1);
}

BOOST_AUTO_TEST_CASE( ZETA_1 )
{
  test_update("/real_update_zeta_input_1", 2);
}

BOOST_AUTO_TEST_SUITE_END()
