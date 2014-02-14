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
#include "../physics/prng.hpp"
#include "../hardware/device.hpp"
#include "../hardware/code/correlator_staggered.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE_CORRELATORS_STAGGERED
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

void fill_sf_with_zero(su3vec * sf_in, int size)
{
  for(int i = 0; i < size; ++i) {
    sf_in[i].e0 = hmc_complex_zero;
    sf_in[i].e1 = hmc_complex_zero;
    sf_in[i].e2 = hmc_complex_zero;
  }
  return;
}

//This function sums the real and imaginary parts of all su3vec contained in sf_in
hmc_float count_sf(su3vec * sf_in, int size)
{
  hmc_float sum = 0.;
  for (int i=0; i<size; i++){
    sum +=
        sf_in[i].e0.re + sf_in[i].e0.im 
      + sf_in[i].e1.re + sf_in[i].e1.im 
      + sf_in[i].e2.re + sf_in[i].e2.im;
  }
  return sum;
}

//The following two function return the sum of the square deviation (frome the mean) of the numbers
//in sf_in. To get the variance, i.e. the mean square deviation, the result
//must be divided by the number of numbers summed (see test_sf_gaussian_staggered in this file).
hmc_float calc_var(hmc_float in, hmc_float mean){
  return (in - mean) * (in - mean);
}

hmc_float calc_var_sf(su3vec * sf_in, int size, hmc_float sum){
  hmc_float var = 0.;
  for(int k=0; k<size; k++){
    var +=
        calc_var(sf_in[k].e0.re, sum) 
      + calc_var(sf_in[k].e0.im, sum) 
      + calc_var(sf_in[k].e1.re, sum)
      + calc_var(sf_in[k].e1.im, sum) 
      + calc_var(sf_in[k].e2.re, sum) 
      + calc_var(sf_in[k].e2.im, sum);
  }
  return var;
}

//This function fills the field sf_in in the following way
// eo==true  ---> sf_in[even]=ONE  and sf_in[odd]=ZERO
// eo==false ---> sf_in[even]=ZERO and sf_in[odd]=ONE
void fill_sf_with_one_eo(su3vec * sf_in, int size, bool eo, meta::Inputparameters &params)
{
  int ns = params.get_nspace();
  int nt = params.get_ntime();
  int x,y,z,t;
  for (x = 0; x<ns; x++){
    for (y = 0; y<ns; y++){
      for (z = 0; z<ns; z++){
	for (t = 0; t<nt; t++){
	  int coord[4];
	  coord[0] = t;
	  coord[1] = x;
	  coord[2] = y;
	  coord[3] = z;
	  int nspace =  get_nspace(coord, params);
	  int global_pos = get_global_pos(nspace, t, params);
	  //This if should be unnecessary if size==ns*ns*ns*nt
	  if (global_pos > size)
	    break;
	  hmc_complex content;
	  if ((x+y+z+t) %2 == 0){
	    if (eo)
	      content = hmc_complex_one;
	    else
	      content = hmc_complex_zero;
	  }
	  else{
	    if (eo)
	      content = hmc_complex_zero;
	    else
	      content = hmc_complex_one;
	  }
	  
	  sf_in[global_pos].e0 = content;
	  sf_in[global_pos].e1 = content;
	  sf_in[global_pos].e2 = content;
	}}}}
  return;
}

//This function works in the following way
// eo==true  ---> sum of all components of sf_in[even] is returned
// eo==false ---> sum of all components of sf_in[odd]  is returned
hmc_float count_sf_eo(su3vec * sf_in, int size, bool eo, meta::Inputparameters &params)
{
  int ns = params.get_nspace();
  int nt = params.get_ntime();
  int x,y,z,t;
  hmc_float sum = 0.;
  for (x = 0; x<ns; x++){
    for (y = 0; y<ns; y++){
      for (z = 0; z<ns; z++){
	for (t = 0; t<nt; t++){
	  int coord[4];
	  coord[0] = t;
	  coord[1] = x;
	  coord[2] = y;
	  coord[3] = z;
	  int nspace =  get_nspace(coord, params);
	  int global_pos = get_global_pos(nspace, t, params);
	  //This if should be unnecessary if size==ns*ns*ns*nt
	  if (global_pos > size)
	    break;
	  if (
	      ( eo ==true && (x+y+z+t) %2 == 0) ||
	      ( eo ==false &&  (x+y+z+t) %2 == 1 )
	      )
	    {
	      int i = global_pos;
	      sum +=
		sf_in[i].e0.re+ sf_in[i].e0.im 
		+sf_in[i].e1.re+ sf_in[i].e1.im 
		+sf_in[i].e2.re+ sf_in[i].e2.im;
	    }
	  else{
	    continue;
	  }
	}}}}
  return sum;
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


/********************************************************************************/

void test_build(std::string inputfile)
{
	logger.info() << "build opencl_module_correlators_staggered";
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	for(auto device: system.get_devices()) {
		device->get_correlator_staggered_code();
	}
	BOOST_MESSAGE("Test done");
}

void test_src_volume(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "create_staggered_eo_volume_source";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);

	physics::PRNG prng(system);
	auto * device = system.get_devices().at(0)->get_correlator_staggered_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_eoprec_spinorfieldsize(params);
	const SU3vec out(NUM_ELEMENTS_SF, device->get_device());

	//The number of times the kernel is run
	int iterations = params.get_integrationsteps(0);
	
	su3vec * sf_out;
	sf_out = new su3vec[NUM_ELEMENTS_SF * iterations];
	BOOST_REQUIRE(sf_out);
	auto prng_buf = prng.get_buffers().at(0);

	hmc_float sum = 0;
	for (int i = 0; i< iterations; i++){
	  if(i%200==0)logger.info() << "Run kernel for the " << i << "th time";
	  out.clear();
	  device->create_volume_source_stagg_eoprec_device(&out, prng_buf);
	  out.dump(&sf_out[i*NUM_ELEMENTS_SF]);
	  //Here we sum the entries to calculate the mean later
	  sum += count_sf(&sf_out[i*NUM_ELEMENTS_SF], NUM_ELEMENTS_SF);
	}
	logger.info() << "result: mean";
	hmc_float cpu_res;
	//sum is the sum of iterations*NUM_ELEMENTS_SF*6 real numbers
	if(params.get_sourcecontent() == meta::Inputparameters::z2)
	  sum = sum/iterations/NUM_ELEMENTS_SF/3;//because immaginary part is not randomly drawn, it is 0.0 always
	else
	  sum = sum/iterations/NUM_ELEMENTS_SF/6;
	cpu_res= sum;
	logger.info() << cpu_res;

	if(params.get_read_multiple_configs()  == false){
	  //Calc std derivation
	  hmc_float var=0.;
	  for (int i=0; i<iterations; i++){
	    var += calc_var_sf(&sf_out[i*NUM_ELEMENTS_SF], NUM_ELEMENTS_SF, sum);
	  }
	  //var is the sum of iterations*NUM_ELEMENTS_SF*6 square deviations
	  if(params.get_sourcecontent() == meta::Inputparameters::z2)
	    var=var/iterations/NUM_ELEMENTS_SF/3;//because immaginary part is not randomly drawn, it is 0.0 always
	  else
	    var=var/iterations/NUM_ELEMENTS_SF/6;
	  
	  cpu_res = sqrt(var);
	  logger.info() << "result: variance";
	  logger.info() << cpu_res;
	}

	if(params.get_sourcecontent() == meta::Inputparameters::z2){
	  if(params.get_read_multiple_configs()  == false)
	    testFloatAgainstInputparameters(cpu_res, params);
	  else
	    testFloatSizeAgainstInputparameters(cpu_res, params);	  
	}else if(params.get_sourcecontent() == meta::Inputparameters::one){
	  testFloatAgainstInputparameters(cpu_res, params);
	} else{
	  testFloatSizeAgainstInputparameters(cpu_res, params);
	}
	BOOST_MESSAGE("Test done");
}



BOOST_AUTO_TEST_SUITE(BUILD)

BOOST_AUTO_TEST_CASE( BUILD_1 )
{
  test_build("/opencl_module_correlators_staggered_build_input_1"); //just to have an input file
}

BOOST_AUTO_TEST_CASE( BUILD_2 )
{
  test_build("/opencl_module_correlators_staggered_build_input_2"); //just to have an input file
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(SRC_VOLUME)

BOOST_AUTO_TEST_CASE( SRC_VOLUME_1 )
{
  test_src_volume("/src_volume_staggered_eo_input_1");
}

BOOST_AUTO_TEST_CASE( SRC_VOLUME_2 )
{
  test_src_volume("/src_volume_staggered_eo_input_2");
}

BOOST_AUTO_TEST_CASE( SRC_VOLUME_3 )
{
  test_src_volume("/src_volume_staggered_eo_input_3");
}

BOOST_AUTO_TEST_CASE( SRC_VOLUME_4 )
{
  test_src_volume("/src_volume_staggered_eo_input_4");
}

BOOST_AUTO_TEST_CASE( SRC_VOLUME_5 )
{
  test_src_volume("/src_volume_staggered_eo_input_5");
}

BOOST_AUTO_TEST_CASE( SRC_VOLUME_6 )
{
  test_src_volume("/src_volume_staggered_eo_input_6");
}

BOOST_AUTO_TEST_CASE( SRC_VOLUME_7 )
{
  test_src_volume("/src_volume_staggered_eo_input_7");
}

BOOST_AUTO_TEST_SUITE_END()

