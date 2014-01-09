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
// #include "../host_random.h"
#include "../physics/prng.hpp"
#include "../hardware/device.hpp"
#include "../hardware/code/real.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE_COMPLEX
#include <boost/test/unit_test.hpp>

//some functionality
#include "test_util.h"

void test_build(std::string inputfile)
{
	logger.info() << "build opencl_module_real";
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	for(auto device: system.get_devices()) {
		device->get_real_code();
	}
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
	meta::Inputparameters params = create_parameters(inputfile);
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
	test_build("/opencl_module_spinors_build_input_1"); //Just to use an inputfile
}

BOOST_AUTO_TEST_CASE( BUILD_2 )
{
	test_build("/opencl_module_spinors_build_input_2"); //Just to use an inputfile
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
