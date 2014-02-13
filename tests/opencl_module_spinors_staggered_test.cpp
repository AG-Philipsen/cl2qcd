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
#include "../hardware/code/spinors_staggered.hpp"
//spinors.hpp needed for get_spinorfieldsize and get_eoprec_spinorfieldsize
#include "../hardware/code/spinors.hpp" 
#include <numeric>      // std::accumulate

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE_SPINORS_STAGGERED
#include <boost/test/unit_test.hpp>

//some functionality
#include "test_util.h"
#include "test_util_staggered.h"
#include "Kolmogorov_Smirnov.h"
#include "Normal_RNG_tests.h"

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
	logger.info() << "build opencl_module_spinors_staggered";
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	for(auto device: system.get_devices()) {
		device->get_spinor_staggered_code();
	}
	BOOST_MESSAGE("Test done");
}

void test_sf_squarenorm_staggered(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "global_squarenorm";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto device = system.get_devices().at(0)->get_spinor_staggered_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_spinorfieldsize(params);	
	const Plain<su3vec> in(NUM_ELEMENTS_SF, device->get_device());
	Plain<hmc_float> sqnorm(1, device->get_device());

	su3vec * sf_in;
	sf_in = new su3vec[NUM_ELEMENTS_SF];
	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	else fill_sf_with_random(sf_in, NUM_ELEMENTS_SF);
	BOOST_REQUIRE(sf_in);

        //The following three lines are to be used to produce the ref_vec file needed to get the ref_value
        //---> Comment them out when the reference values have been obtained! 
        /*
        print_staggeredfield_to_textfile("ref_vec_sq",sf_in,params); 
        logger.info() << "Produced the ref_vec text file with the staggered field for the ref. code. Returning...";   
        return;
	// */

	in.load(sf_in);

	logger.info() << "Run kernel";
	logger.info() << "result:";
	hmc_float cpu_res;
	device->set_float_to_global_squarenorm_device(&in, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_sf_scalar_product_staggered(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "scalar_product";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto device = system.get_devices().at(0)->get_spinor_staggered_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_spinorfieldsize(params);
	const Plain<su3vec> in(NUM_ELEMENTS_SF, device->get_device());
	const Plain<su3vec> in2(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_complex> sqnorm(1, device->get_device());

	su3vec * sf_in;
	sf_in = new su3vec[NUM_ELEMENTS_SF];
	su3vec * sf_in2;
	sf_in2 = new su3vec[NUM_ELEMENTS_SF];
	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) {
	  fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	  fill_sf_with_one(sf_in2, NUM_ELEMENTS_SF);
	}
	else {
	  fill_sf_with_random(sf_in, NUM_ELEMENTS_SF, 123);
	  fill_sf_with_random(sf_in2, NUM_ELEMENTS_SF, 456);
	}
	BOOST_REQUIRE(sf_in);
	BOOST_REQUIRE(sf_in2);

        //The following five lines are to be used to produce the ref_vec file needed to get the ref_value
        //---> Comment them out when the reference values have been obtained! 
        /*
        print_staggeredfield_to_textfile("ref_vec_sp1",sf_in,params); 
        logger.info() << "Produced the ref_vec_sp1 text file with the staggered field for the ref. code.";   
        print_staggeredfield_to_textfile("ref_vec_sp2",sf_in2,params); 
        logger.info() << "Produced the ref_vec_sp2 text file with the staggered field for the ref. code. Returning...";   
        return;
	// */

	in.load(sf_in);
	in2.load(sf_in2);

	logger.info() << "Run kernel";
	logger.info() << "result:";
	hmc_complex cpu_res_tmp;
	device->set_complex_to_scalar_product_device(&in, &in2, &sqnorm);
	sqnorm.dump(&cpu_res_tmp);
	hmc_float cpu_res = cpu_res_tmp.re + cpu_res_tmp.im;
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_sf_cold(std::string inputfile, bool switcher)
{
	//switcher decides if the sf is set to cold or zero   
	using namespace hardware::buffers;

	std::string kernelName;
	if(switcher)
	  kernelName = "set_cold_spinorfield_stagg";
	else
	  kernelName = "set_zero_spinorfield_stagg";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto device = system.get_devices().at(0)->get_spinor_staggered_code();
	
	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_spinorfieldsize(params);
	const Plain<su3vec> in(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	logger.info() << "Run kernel";
        if(switcher)
	  device->set_cold_spinorfield_device(&in);
        else
          device->set_zero_spinorfield_device(&in);

	logger.info() << "result:";
	hmc_float cpu_res;
	device->set_float_to_global_squarenorm_device(&in, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_sf_sax_staggered(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "sax_staggered";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto device = system.get_devices().at(0)->get_spinor_staggered_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_spinorfieldsize(params);
	const Plain<su3vec> in(NUM_ELEMENTS_SF, device->get_device());
	const Plain<su3vec> out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	hardware::buffers::Plain<hmc_complex> alpha(1, device->get_device());

	hmc_complex alpha_host = {params.get_beta(), params.get_rho()};
	logger.info() << "Use alpha = (" << alpha_host.re << ","<< alpha_host.im <<")";

	su3vec * sf_in;
	sf_in = new su3vec[NUM_ELEMENTS_SF];
	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg)
	  fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	else
	  fill_sf_with_random(sf_in, NUM_ELEMENTS_SF, 123);
	BOOST_REQUIRE(sf_in);

	//The following three lines are to be used to produce the ref_vec file needed to get the ref_value
        //---> Comment them out when the reference values have been obtained! 
        /*
        print_staggeredfield_to_textfile("ref_vec_sax",sf_in,params); 
        logger.info() << "Produced the ref_vec text file with the staggered field for the ref. code. Returning...";   
        return;
	// */
	
	in.load(sf_in);
	alpha.load(&alpha_host);

	logger.info() << "Run kernel";
	device->sax_device(&in, &alpha, &out);

	logger.info() << "result:";
	hmc_float cpu_res;
	device->set_float_to_global_squarenorm_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_sf_saxpy_staggered(std::string inputfile, bool switcher)
{
  //switcher chooses between saxpy_stagg and saxpy_stagg_arg kernel, which have the same functionality
	using namespace hardware::buffers;

	std::string kernelName;
	if( switcher)
	  kernelName = "saxpy_staggered";
	else {
	  //kernelName = "saxpy_staggered_arg";
	  ;
	}	  
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto device = system.get_devices().at(0)->get_spinor_staggered_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_spinorfieldsize(params);
	const Plain<su3vec> in(NUM_ELEMENTS_SF, device->get_device());
	const Plain<su3vec> in2(NUM_ELEMENTS_SF, device->get_device());
	const Plain<su3vec> out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	hardware::buffers::Plain<hmc_complex> alpha(1, device->get_device());

	hmc_complex alpha_host = {params.get_beta(), params.get_rho()};
	logger.info() << "Use alpha = (" << alpha_host.re << ","<< alpha_host.im <<")";

	su3vec * sf_in;
	su3vec * sf_in2;
	sf_in = new su3vec[NUM_ELEMENTS_SF];
	sf_in2 = new su3vec[NUM_ELEMENTS_SF];
	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) {
	  fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	  fill_sf_with_one(sf_in2, NUM_ELEMENTS_SF);
	}
	else {
	  fill_sf_with_random(sf_in, NUM_ELEMENTS_SF, 123);
	  fill_sf_with_random(sf_in2, NUM_ELEMENTS_SF, 456);
	}
	BOOST_REQUIRE(sf_in);
	BOOST_REQUIRE(sf_in2);

	//The following five lines are to be used to produce the ref_vec file needed to get the ref_value
        //---> Comment them out when the reference values have been obtained! 
        /*
        print_staggeredfield_to_textfile("ref_vec_saxpy1",sf_in,params); 
        logger.info() << "Produced the ref_vec_saxpy1 text file with the staggered field for the ref. code.";   
        print_staggeredfield_to_textfile("ref_vec_saxpy2",sf_in2,params); 
        logger.info() << "Produced the ref_vec_saxpy2 text file with the staggered field for the ref. code. Returning...";   
        return;
	// */
	
	in.load(sf_in);
	in2.load(sf_in2);
	alpha.load(&alpha_host);

	logger.info() << "Run kernel";
	if (switcher)
	  device->saxpy_device(&in, &in2, &alpha, &out);
	else 
	  logger.error() << "The kernel saxpy_stagg_arg doesn't exist yet, this is a meaningless test!";
	  //device->saxpy_device(&in, &in2, alpha_host, &out);

	logger.info() << "result:";
	hmc_float cpu_res;
	device->set_float_to_global_squarenorm_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_sf_saxpbypz_staggered(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "saxpbypz_stagg";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto device = system.get_devices().at(0)->get_spinor_staggered_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_spinorfieldsize(params);
	const Plain<su3vec> in(NUM_ELEMENTS_SF, device->get_device());
	const Plain<su3vec> in2(NUM_ELEMENTS_SF, device->get_device());
	const Plain<su3vec> in3(NUM_ELEMENTS_SF, device->get_device());
	const Plain<su3vec> out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	hardware::buffers::Plain<hmc_complex> alpha(1, device->get_device());
	hardware::buffers::Plain<hmc_complex> beta(1, device->get_device());

	hmc_complex alpha_host = {params.get_beta(), params.get_rho()};
	hmc_complex beta_host = {params.get_kappa(), params.get_mu()};
	logger.info() << "Use alpha = (" << alpha_host.re << ","<< alpha_host.im <<")";
	logger.info() << "Use beta = (" << beta_host.re << ","<< beta_host.im <<")";

	su3vec * sf_in;
	su3vec * sf_in2;
	su3vec * sf_in3;
	sf_in = new su3vec[NUM_ELEMENTS_SF];
	sf_in2 = new su3vec[NUM_ELEMENTS_SF];
	sf_in3 = new su3vec[NUM_ELEMENTS_SF];
	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) {
	  fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	  fill_sf_with_one(sf_in2, NUM_ELEMENTS_SF);
	  fill_sf_with_one(sf_in3, NUM_ELEMENTS_SF);
	}
	else {
	  fill_sf_with_random(sf_in, NUM_ELEMENTS_SF, 123);
	  fill_sf_with_random(sf_in2, NUM_ELEMENTS_SF, 456);
	  fill_sf_with_random(sf_in3, NUM_ELEMENTS_SF, 789);
	}
	BOOST_REQUIRE(sf_in);
	BOOST_REQUIRE(sf_in2);
	BOOST_REQUIRE(sf_in3);

	//The following seven lines are to be used to produce the ref_vec file needed to get the ref_value
        //---> Comment them out when the reference values have been obtained! 
	/*
        print_staggeredfield_to_textfile("ref_vec_saxpbypz1",sf_in,params); 
        logger.info() << "Produced the ref_vec_saxpbypz1 text file with the staggered field for the ref. code."; 
	print_staggeredfield_to_textfile("ref_vec_saxpbypz2",sf_in2,params); 
        logger.info() << "Produced the ref_vec_saxpbypz2 text file with the staggered field for the ref. code.";  
        print_staggeredfield_to_textfile("ref_vec_saxpbypz3",sf_in3,params); 
        logger.info() << "Produced the ref_vec_saxpbypz3 text file with the staggered field for the ref. code. Returning...";   
        return;
	// */
	
	in.load(sf_in);
	in2.load(sf_in2);
	in3.load(sf_in3);
	alpha.load(&alpha_host);
	beta.load(&beta_host);

	logger.info() << "Run kernel";
	device->saxpbypz_device(&in, &in2, &in3, &alpha, &beta, &out);

	logger.info() << "result:";
	hmc_float cpu_res;
	device->set_float_to_global_squarenorm_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_sf_gaussian_staggered(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "set_gaussian_spinorfield_stagg";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);

	physics::PRNG prng(system);
	auto device = system.get_devices().at(0)->get_spinor_staggered_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_spinorfieldsize(params);
	const Plain<su3vec> out(NUM_ELEMENTS_SF, device->get_device());

	int iterations = params.get_integrationsteps(0);

	su3vec * sf_out;
	sf_out = new su3vec[NUM_ELEMENTS_SF * iterations];
	BOOST_REQUIRE(sf_out);

	auto prng_buf = prng.get_buffers().at(0);

	hmc_float sum = 0;
	for (int i = 0; i< iterations; i++){
	  if(i%100==0)logger.info() << "Run kernel for the " << i << "th time";
	  device->set_gaussian_spinorfield_device(&out, prng_buf);
	  out.dump(&sf_out[i*NUM_ELEMENTS_SF]);
	  //Here we sum the entries to calculate the mean later
	  sum += count_sf(&sf_out[i*NUM_ELEMENTS_SF], NUM_ELEMENTS_SF);
	}
	
	logger.info() << "result: mean";
	hmc_float cpu_res = 0.;
	//sum is the sum of iterations*NUM_ELEMENTS_SF*6 real gaussian numbers
	sum = sum/iterations/NUM_ELEMENTS_SF/6;	
	cpu_res= sum;
	logger.info() << cpu_res;
	
	if(params.get_read_multiple_configs()  == false){
	  //CP: calc std derivation
	  hmc_float var=0.;
	  for (int i=0; i<iterations; i++){
	    var += calc_var_sf(&sf_out[i*NUM_ELEMENTS_SF], NUM_ELEMENTS_SF, sum);
	  }
	  //var is the sum of iterations*NUM_ELEMENTS_SF*6 square deviations
	  var=var/iterations/NUM_ELEMENTS_SF/6;
	  
	  cpu_res = sqrt(var);
	  logger.info() << "result: variance";
	  logger.info() << cpu_res;
	}
	
	//The Kolmogorov_Smirnov test requires set of n samples with n around 1000
	//(to big n and to small n are not good choices for this test)
	//So we can consider each kernel result as a set (NUM_ELEMENTS_SF*6=1536 for a 4^4 lattice)
	vector<vector<hmc_float>> samples;
	vector<hmc_float> tmp;
	vector<hmc_float> tmp2;
	for(int i=0; i<iterations; i++){
	  vector<hmc_float> tmp;
	  for(uint j=0; j<NUM_ELEMENTS_SF; j++){
	    tmp2=reals_from_su3vec(sf_out[i*NUM_ELEMENTS_SF+j]);
	    tmp.insert(tmp.end(),tmp2.begin(),tmp2.end());
	    tmp2.clear();
	  }
	  samples.push_back(tmp);
	  tmp.clear();
	}
	logger.info() << "Running Kolmogorov_Smirnov test (it should take half a minute)...";
	logger.info() << "Kolmogorov_Smirnov frequency (of K+): " << std::setprecision(16) << Kolmogorov_Smirnov(samples,0.,sqrt(0.5)) << " ---> It should be close 0.98";
	
	if(params.get_read_multiple_configs()==true){
	  //Let us use the same sets of samples to make the mean test to 1,2 and 3 sigma
	  //Note that in the test BOOST_CHECK is used.
	  mean_test_multiple_set(samples,2.,0.,sqrt(0.5));
	  mean_test_multiple_set(samples,3.,0.,sqrt(0.5));
	  mean_test_multiple_set(samples,4.,0.,sqrt(0.5));
	}else{
	  //Let us use the same sets of samples to make the variance test to 1,2 and 3 sigma
	  //Note that in the test BOOST_CHECK is used.
	  variance_test_multiple_set(samples,2.,sqrt(0.5));
	  variance_test_multiple_set(samples,3.,sqrt(0.5));
	  variance_test_multiple_set(samples,4.,sqrt(0.5));
	}
	
	//Here we test if cpu_res is smaller than ref_value: in this case the test passes
	testFloatSizeAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////           EVEN-ODD PRECONDITIONING             //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void test_sf_convert_to_eo_staggered(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "convert_to_eoprec_staggered";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto device = system.get_devices().at(0)->get_spinor_staggered_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF_EO = hardware::code::get_eoprec_spinorfieldsize(params);
	size_t NUM_ELEMENTS_SF = hardware::code::get_spinorfieldsize(params);
	const Plain<su3vec> in(NUM_ELEMENTS_SF, device->get_device());
	const SU3vec in2(NUM_ELEMENTS_SF_EO, device->get_device());
	const SU3vec in3(NUM_ELEMENTS_SF_EO, device->get_device());
	const Plain<su3vec> out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	su3vec * sf_in;
	sf_in = new su3vec[NUM_ELEMENTS_SF];
	if(params.get_read_multiple_configs() )
	  fill_sf_with_one_eo(sf_in, NUM_ELEMENTS_SF, true, params);
	else
	  fill_sf_with_one_eo(sf_in, NUM_ELEMENTS_SF, false, params);
	BOOST_REQUIRE(sf_in);

	in.load(sf_in);
	
	logger.info() << "Run kernel";
	device->convert_to_eoprec_device(&in2, &in3, &in);

	logger.info() << "result:";
	hmc_float cpu_res;
	if(params.get_read_multiple_configs() ){
	  device->set_float_to_global_squarenorm_eoprec_device(&in3, &sqnorm);
	  sqnorm.dump(&cpu_res);
	  logger.info() << cpu_res;
	  //CP: this must be zero since only the even sites should be filled!
	  BOOST_REQUIRE_CLOSE(cpu_res, 0., params.get_solver_prec());
	  device->set_float_to_global_squarenorm_eoprec_device(&in2, &sqnorm);
	  sqnorm.dump(&cpu_res);
	  logger.info() << cpu_res;
	}
	else{
	  device->set_float_to_global_squarenorm_eoprec_device(&in2, &sqnorm);
	  sqnorm.dump(&cpu_res);
	  logger.info() << cpu_res;
	  //CP: this must be zero since only the odd sites should be filled!
	  BOOST_REQUIRE_CLOSE(cpu_res, 0., params.get_solver_prec());
	  device->set_float_to_global_squarenorm_eoprec_device(&in3, &sqnorm);
	  sqnorm.dump(&cpu_res);
	  logger.info() << cpu_res;
	}

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_sf_convert_from_eo_staggered(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "convert_from_eoprec_staggered";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto device = system.get_devices().at(0)->get_spinor_staggered_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF_EO = hardware::code::get_eoprec_spinorfieldsize(params);
	size_t NUM_ELEMENTS_SF = hardware::code::get_spinorfieldsize(params);
	const SU3vec in2(NUM_ELEMENTS_SF_EO, device->get_device());
	const SU3vec in3(NUM_ELEMENTS_SF_EO, device->get_device());
	const Plain<su3vec> out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	su3vec * sf_eo1;
	su3vec * sf_eo2;
	su3vec * sf_out;
	sf_eo1 = new su3vec[NUM_ELEMENTS_SF_EO];
	sf_eo2 = new su3vec[NUM_ELEMENTS_SF_EO];
	sf_out = new su3vec[NUM_ELEMENTS_SF];
	if(params.get_read_multiple_configs() ){
	  fill_sf_with_one(sf_eo1, NUM_ELEMENTS_SF_EO);
	  fill_sf_with_zero(sf_eo2, NUM_ELEMENTS_SF_EO);
	}else{
	  fill_sf_with_zero(sf_eo1, NUM_ELEMENTS_SF_EO);
	  fill_sf_with_one(sf_eo2, NUM_ELEMENTS_SF_EO);
	}
	BOOST_REQUIRE(sf_eo1);
	BOOST_REQUIRE(sf_eo2);

	in2.load(sf_eo1);
	in3.load(sf_eo2);

	logger.info() << "Run kernel";
	device->convert_from_eoprec_device(&in2, &in3, &out);

	out.dump(sf_out);

	logger.info() << "result:";
	hmc_float cpu_res;
	if(params.get_read_multiple_configs() ){
	  cpu_res= count_sf_eo(sf_out, NUM_ELEMENTS_SF, true, params);	
	  logger.info() << cpu_res;
	}
	else{
	  cpu_res= count_sf_eo(sf_out, NUM_ELEMENTS_SF, false, params);	
	  logger.info() << cpu_res;
	}

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}


void test_sf_squarenorm_staggered_eo(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "global_squarenorm_staggered_eoprec";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto device = system.get_devices().at(0)->get_spinor_staggered_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_eoprec_spinorfieldsize(params);
	const SU3vec in(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	su3vec * sf_in;
	sf_in = new su3vec[NUM_ELEMENTS_SF];
	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	else fill_sf_with_random(sf_in, NUM_ELEMENTS_SF);
	BOOST_REQUIRE(sf_in);
	
	//The following three lines are to be used to produce the ref_vec file needed to get the ref_value
        //---> Comment them out when the reference values have been obtained! 
        /*
        print_staggeredfield_eo_to_textfile("ref_vec_sq_eo",sf_in,params); 
        logger.info() << "Produced the ref_vec_sq_eo text file with the staggered field for the ref. code. Returning...";   
        return;
	// */
	
	in.load(sf_in);
	
	logger.info() << "Run kernel";
	logger.info() << "result:";
	hmc_float cpu_res;
	device->set_float_to_global_squarenorm_eoprec_device(&in, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_sf_cold_staggered_eo(std::string inputfile, bool switcher)
{
  //switcher decides if the sf is set to cold or zero
	using namespace hardware::buffers;

	std::string kernelName;
	if(switcher)
	  kernelName = "set_cold_spinorfield_stagg_eoprec";
	else
	  kernelName = "set_zero_spinorfield_stagg_eoprec";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto device = system.get_devices().at(0)->get_spinor_staggered_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_eoprec_spinorfieldsize(params);
	const SU3vec in(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	logger.info() << "Run kernel";
	if(switcher)
	  device->set_cold_spinorfield_eoprec_device(&in);
	else
	  device->set_zero_spinorfield_eoprec_device(&in);
	logger.info() << "result:";
	hmc_float cpu_res;
	device->set_float_to_global_squarenorm_eoprec_device(&in, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_sf_scalar_product_staggered_eo(std::string inputfile, bool real_part=false)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "scalar_product_eoprec_staggered";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto device = system.get_devices().at(0)->get_spinor_staggered_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_eoprec_spinorfieldsize(params);
	const SU3vec in(NUM_ELEMENTS_SF, device->get_device());
	const SU3vec in2(NUM_ELEMENTS_SF, device->get_device());
	//Here I waste a bit of memory but in a test this is not so serious
	hardware::buffers::Plain<hmc_complex> sqnorm(1, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm_real(1, device->get_device());

	su3vec * sf_in;
	sf_in = new su3vec[NUM_ELEMENTS_SF];
	su3vec * sf_in2;
	sf_in2 = new su3vec[NUM_ELEMENTS_SF];
	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) {
	  fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	  fill_sf_with_one(sf_in2, NUM_ELEMENTS_SF);
	}
	else {
	  fill_sf_with_random(sf_in, NUM_ELEMENTS_SF, 123);
	  fill_sf_with_random(sf_in2, NUM_ELEMENTS_SF, 456);
	}
	BOOST_REQUIRE(sf_in);
	BOOST_REQUIRE(sf_in2);

	//The following five lines are to be used to produce the ref_vec file needed to get the ref_value
        //---> Comment them out when the reference values have been obtained! 
        /*
        print_staggeredfield_eo_to_textfile("ref_vec_sp1",sf_in,params); 
        logger.info() << "Produced the ref_vec_sp1 text file with the staggered field for the ref. code.";   
        print_staggeredfield_eo_to_textfile("ref_vec_sp2",sf_in2,params); 
        logger.info() << "Produced the ref_vec_sp2 text file with the staggered field for the ref. code. Returning...";   
        return;
	// */
	
	in.load(sf_in);
	in2.load(sf_in2);

	logger.info() << "Run kernel";
	logger.info() << "result:";
	hmc_complex cpu_res_tmp;
	hmc_float cpu_res;
	if(real_part == false){
	  device->set_complex_to_scalar_product_eoprec_device(&in, &in2, &sqnorm);
	  sqnorm.dump(&cpu_res_tmp);
	  cpu_res = cpu_res_tmp.re + cpu_res_tmp.im;
	}else{
	  device->set_float_to_scalar_product_real_part_eoprec_device(&in, &in2, &sqnorm_real);
	  sqnorm_real.dump(&cpu_res);
	}
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_sf_sax_staggered_eo(std::string inputfile, int switcher=0)
{
 //switcher chooses between sax and sax_arg kernel (real or complex), which have the same functionality
	using namespace hardware::buffers;

	std::string kernelName;
	if(switcher==0)
	  kernelName = "sax_cplx_staggered_eoprec";
	if(switcher==1)
	  kernelName = "sax_cplx_arg_staggered_eoprec";
	if(switcher==2)
	  kernelName = "sax_real_staggered_eoprec";
	if(switcher==3)
	  kernelName = "sax_real_arg_staggered_eoprec";
	if(switcher==4)
	  kernelName = "sax_real_vec_staggered_eoprec";
	
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto device = system.get_devices().at(0)->get_spinor_staggered_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_eoprec_spinorfieldsize(params);
	const SU3vec in(NUM_ELEMENTS_SF, device->get_device());
	const SU3vec out(NUM_ELEMENTS_SF, device->get_device());
	
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	//Here we waste a bit of memory, but in the test is not a problem!!
	//Used if switcher==4
	hardware::buffers::Plain<hmc_float> alpha_real_vec(5, device->get_device());
	std::vector<hmc_float> alpha_host_real_vec(5, params.get_beta());
	const int index_alpha = 3;
	//Used if switcher==2 || ==3
	hardware::buffers::Plain<hmc_float> alpha_real(1, device->get_device());
	hmc_float alpha_host_real = params.get_beta();
	//Used if switcher==0 || ==1
	hardware::buffers::Plain<hmc_complex> alpha(1, device->get_device());
	hmc_complex alpha_host = {params.get_beta(), params.get_rho()};

	su3vec * sf_in;
	sf_in = new su3vec[NUM_ELEMENTS_SF];
	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) {
	  fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	}
	else {
	  fill_sf_with_random(sf_in, NUM_ELEMENTS_SF, 123);
	}
	BOOST_REQUIRE(sf_in);

	//The following three lines are to be used to produce the ref_vec file needed to get the ref_value
        //---> Comment them out when the reference values have been obtained! 
        /*
        print_staggeredfield_eo_to_textfile("ref_vec_sax_eo",sf_in,params); 
        logger.info() << "Produced the ref_vec text file with the staggered field for the ref. code. Returning...";   
        return;
	// */
	
	in.load(sf_in);
	if(switcher==4){
	  logger.info() << "Use alpha[" << index_alpha << "] = " << alpha_host_real_vec[index_alpha];
	  alpha_real_vec.load(&alpha_host_real_vec[0]); 
	}else if(switcher==2 || switcher==3){
	  logger.info() << "Use alpha = " << alpha_host_real;
	  alpha_real.load(&alpha_host_real);
	}else{
	  logger.info() << "Use alpha = (" << alpha_host.re << ","<< alpha_host.im <<")";
	  alpha.load(&alpha_host);
	}
	

	logger.info() << "Run kernel";
	if(switcher==0)
	  device->sax_eoprec_device(&in, &alpha, &out);
	if(switcher==1)
	  device->sax_eoprec_device(&in, alpha_host, &out);
	if(switcher==2)
	  device->sax_eoprec_device(&in, &alpha_real, &out);
	if(switcher==3)
	  device->sax_eoprec_device(&in, alpha_host_real, &out);
	if(switcher==4)
	  device->sax_eoprec_device(&in, &alpha_real_vec, index_alpha, &out);
	
	logger.info() << "result:";
	hmc_float cpu_res;
	device->set_float_to_global_squarenorm_eoprec_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_sf_saxpy_staggered_eo(std::string inputfile, int switcher=0)
{
  //switcher chooses between saxpy and saxpy_arg kernel (real or complex), which have the same functionality
	using namespace hardware::buffers;

	std::string kernelName;
	if(switcher==0)
	  kernelName = "saxpy_cplx_staggered_eoprec";
	if(switcher==1)
	  kernelName = "saxpy_cplx_arg_staggered_eoprec";
	if(switcher==2)
	  kernelName = "saxpy_real_staggered_eoprec";
	if(switcher==3)
	  kernelName = "saxpy_real_arg_staggered_eoprec";
	if(switcher==4)
	  kernelName = "saxpy_real_vec_staggered_eoprec";
	
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto device = system.get_devices().at(0)->get_spinor_staggered_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_eoprec_spinorfieldsize(params);
	const SU3vec in(NUM_ELEMENTS_SF, device->get_device());
	const SU3vec in2(NUM_ELEMENTS_SF, device->get_device());
	const SU3vec out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	
	//Here we waste a bit of memory, but in the test is not a problem!!
	//Used if switcher==4
	hardware::buffers::Plain<hmc_float> alpha_real_vec(5, device->get_device());
	std::vector<hmc_float> alpha_host_real_vec(5, params.get_beta());
	const int index_alpha = 3;
	//Used if switcher==2 || ==3
	hardware::buffers::Plain<hmc_float> alpha_real(1, device->get_device());
	hmc_float alpha_host_real = params.get_beta();
	//Used if switcher==0 || ==1
	hardware::buffers::Plain<hmc_complex> alpha(1, device->get_device());
	hmc_complex alpha_host = {params.get_beta(), params.get_rho()};
	
	su3vec * sf_in;
	su3vec * sf_in2;
	sf_in = new su3vec[NUM_ELEMENTS_SF];
	sf_in2 = new su3vec[NUM_ELEMENTS_SF];
	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) {
	  fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	  fill_sf_with_one(sf_in2, NUM_ELEMENTS_SF);
	}
	else {
	  fill_sf_with_random(sf_in, NUM_ELEMENTS_SF, 123);
	  fill_sf_with_random(sf_in2, NUM_ELEMENTS_SF, 456);
	}
	BOOST_REQUIRE(sf_in);
	BOOST_REQUIRE(sf_in2);

	//The following five lines are to be used to produce the ref_vec file needed to get the ref_value
        //---> Comment them out when the reference values have been obtained! 
        /*
        print_staggeredfield_eo_to_textfile("ref_vec_saxpy1_eo",sf_in,params); 
        logger.info() << "Produced the ref_vec_saxpy1_eo text file with the staggered field for the ref. code.";   
        print_staggeredfield_eo_to_textfile("ref_vec_saxpy2_eo",sf_in2,params); 
        logger.info() << "Produced the ref_vec_saxpy2_eo text file with the staggered field for the ref. code. Returning...";   
        return;
	// */
	
	in.load(sf_in);
	in2.load(sf_in2);
	if(switcher==4){
	  logger.info() << "Use alpha[" << index_alpha << "] = " << alpha_host_real_vec[index_alpha];
	  alpha_real_vec.load(&alpha_host_real_vec[0]); 
	}else if(switcher==2 || switcher==3){
	  logger.info() << "Use alpha = " << alpha_host_real;
	  alpha_real.load(&alpha_host_real);
	}else{
	  logger.info() << "Use alpha = (" << alpha_host.re << ","<< alpha_host.im <<")";
	  alpha.load(&alpha_host);
	}
	
	logger.info() << "Run kernel";
	if(switcher==0)
	  device->saxpy_eoprec_device(&in, &in2, &alpha, &out);
	if(switcher==1)
	  device->saxpy_eoprec_device(&in, &in2, alpha_host, &out);
	if(switcher==2)
	  device->saxpy_eoprec_device(&in, &in2, &alpha_real, &out);
	if(switcher==3)
	  device->saxpy_eoprec_device(&in, &in2, alpha_host_real, &out);
	if(switcher==4)
	  device->saxpy_eoprec_device(&in, &in2, &alpha_real_vec, &out);

	logger.info() << "result:";
	hmc_float cpu_res;
	device->set_float_to_global_squarenorm_eoprec_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_sf_saxpby_staggered_eo(std::string inputfile, int switcher=0)
{
	using namespace hardware::buffers;

	std::string kernelName;
	if(switcher==0)
	  kernelName = "saxpby_cplx_staggered_eoprec";
	if(switcher==1)
	  kernelName = "saxpby_cplx_arg_staggered_eoprec";
	if(switcher==2)
	  kernelName = "saxpby_real_staggered_eoprec";
	if(switcher==3)
	  kernelName = "saxpby_real_arg_staggered_eoprec";
	if(switcher==4)
	  kernelName = "saxpby_real_vec_staggered_eoprec";
	
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto device = system.get_devices().at(0)->get_spinor_staggered_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_eoprec_spinorfieldsize(params);
	const SU3vec in(NUM_ELEMENTS_SF, device->get_device());
	const SU3vec in2(NUM_ELEMENTS_SF, device->get_device());
	const SU3vec out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	
	//Here we waste a bit of memory, but in the test is not a problem!!
	//Used if switcher==4
	hardware::buffers::Plain<hmc_float> alpha_real_vec(5, device->get_device());
	std::vector<hmc_float> alpha_host_real_vec(5, params.get_beta());
	const int index_alpha = 3;
	hardware::buffers::Plain<hmc_float> beta_real_vec(5, device->get_device());
	std::vector<hmc_float> beta_host_real_vec(5, params.get_kappa());
	const int index_beta = 2;
	//Used if switcher==2 || ==3
	hardware::buffers::Plain<hmc_float> alpha_real(1, device->get_device());
	hardware::buffers::Plain<hmc_float> beta_real(1, device->get_device());
	hmc_float alpha_host_real = params.get_beta();
	hmc_float beta_host_real = params.get_kappa();
	//Used if switcher==0 || ==1
	hardware::buffers::Plain<hmc_complex> alpha(1, device->get_device());
	hardware::buffers::Plain<hmc_complex> beta(1, device->get_device());
	hmc_complex alpha_host = {params.get_beta(), params.get_rho()};
	hmc_complex beta_host = {params.get_kappa(), params.get_mu()};

	su3vec * sf_in;
	su3vec * sf_in2;
	sf_in = new su3vec[NUM_ELEMENTS_SF];
	sf_in2 = new su3vec[NUM_ELEMENTS_SF];
	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) {
	  fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	  fill_sf_with_one(sf_in2, NUM_ELEMENTS_SF);
	}
	else {
	  fill_sf_with_random(sf_in, NUM_ELEMENTS_SF, 123);
	  fill_sf_with_random(sf_in2, NUM_ELEMENTS_SF, 456);
	}
	BOOST_REQUIRE(sf_in);
	BOOST_REQUIRE(sf_in2);

	//The following seven lines are to be used to produce the ref_vec file needed to get the ref_value
        //---> Comment them out when the reference values have been obtained! 
	/*
        print_staggeredfield_eo_to_textfile("ref_vec_saxpby1_eo",sf_in,params); 
        logger.info() << "Produced the ref_vec_saxpby1_eo text file with the staggered field for the ref. code."; 
	print_staggeredfield_eo_to_textfile("ref_vec_saxpby2_eo",sf_in2,params); 
        logger.info() << "Produced the ref_vec_saxpbypz2_eo text file with the staggered field for the ref. code. Returning...";   
        return;
	// */
	
	in.load(sf_in);
	in2.load(sf_in2);
	
	if(switcher==4){
	  logger.info() << "Use alpha[" << index_alpha << "] = " << alpha_host_real_vec[index_alpha];
	  alpha_real_vec.load(&alpha_host_real_vec[0]); 
	  logger.info() << "Use beta[" << index_beta << "] = " << beta_host_real_vec[index_beta];
	  beta_real_vec.load(&beta_host_real_vec[0]); 
	}else if(switcher==2 || switcher==3){
	  logger.info() << "Use alpha = " << alpha_host_real;
	  logger.info() << "Use beta = " << beta_host_real;
	  alpha_real.load(&alpha_host_real);
	  beta_real.load(&beta_host_real);
	}else{
	  logger.info() << "Use alpha = (" << alpha_host.re << ","<< alpha_host.im <<")";
	  logger.info() << "Use beta = (" << beta_host.re << ","<< beta_host.im <<")";
	  alpha.load(&alpha_host);
	  beta.load(&beta_host);
	}

	logger.info() << "Run kernel";
	if(switcher==0)
	  device->saxpby_eoprec_device(&in, &in2, &alpha, &beta, &out);
	if(switcher==1)
	  device->saxpby_eoprec_device(&in, &in2, alpha_host, beta_host, &out);
	if(switcher==2)
	  device->saxpby_eoprec_device(&in, &in2, &alpha_real, &beta_real, &out);
	if(switcher==3)
	  device->saxpby_eoprec_device(&in, &in2, alpha_host_real, beta_host_real, &out);
	if(switcher==4)
	  device->saxpby_eoprec_device(&in, &in2, &alpha_real_vec, &beta_real_vec, index_alpha, index_beta, &out);
	
	logger.info() << "result:";
	hmc_float cpu_res;
	device->set_float_to_global_squarenorm_eoprec_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_sf_saxpbypz_staggered_eo(std::string inputfile, bool switcher=true)
{
	using namespace hardware::buffers;

	std::string kernelName;
	if(switcher)
	  kernelName = "saxpbypz_cplx_staggered_eoprec";
	else
	  kernelName = "saxpbypz_cplx_arg_staggered_eoprec";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto device = system.get_devices().at(0)->get_spinor_staggered_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_eoprec_spinorfieldsize(params);
	const SU3vec in(NUM_ELEMENTS_SF, device->get_device());
	const SU3vec in2(NUM_ELEMENTS_SF, device->get_device());
	const SU3vec in3(NUM_ELEMENTS_SF, device->get_device());
	const SU3vec out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	hardware::buffers::Plain<hmc_complex> alpha(1, device->get_device());
	hardware::buffers::Plain<hmc_complex> beta(1, device->get_device());

	hmc_complex alpha_host = {params.get_beta(), params.get_rho()};
	hmc_complex beta_host = {params.get_kappa(), params.get_mu()};
	logger.info() << "Use alpha = (" << alpha_host.re << ","<< alpha_host.im <<")";
	logger.info() << "Use beta = (" << beta_host.re << ","<< beta_host.im <<")";

	su3vec * sf_in;
	su3vec * sf_in2;
	su3vec * sf_in3;
	sf_in = new su3vec[NUM_ELEMENTS_SF];
	sf_in2 = new su3vec[NUM_ELEMENTS_SF];
	sf_in3 = new su3vec[NUM_ELEMENTS_SF];
	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) {
	  fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	  fill_sf_with_one(sf_in2, NUM_ELEMENTS_SF);
	  fill_sf_with_one(sf_in3, NUM_ELEMENTS_SF);
	}
	else {
	  fill_sf_with_random(sf_in, NUM_ELEMENTS_SF, 123);
	  fill_sf_with_random(sf_in2, NUM_ELEMENTS_SF, 456);
	  fill_sf_with_random(sf_in3, NUM_ELEMENTS_SF, 789);
	}
	BOOST_REQUIRE(sf_in);
	BOOST_REQUIRE(sf_in2);
	BOOST_REQUIRE(sf_in3);

	//The following seven lines are to be used to produce the ref_vec file needed to get the ref_value
        //---> Comment them out when the reference values have been obtained! 
	/*
        print_staggeredfield_eo_to_textfile("ref_vec_saxpbypz1_eo",sf_in,params); 
        logger.info() << "Produced the ref_vec_saxpbypz1_eo text file with the staggered field for the ref. code."; 
	print_staggeredfield_eo_to_textfile("ref_vec_saxpbypz2_eo",sf_in2,params); 
        logger.info() << "Produced the ref_vec_saxpbypz2_eo text file with the staggered field for the ref. code.";  
        print_staggeredfield_eo_to_textfile("ref_vec_saxpbypz3_eo",sf_in3,params); 
        logger.info() << "Produced the ref_vec_saxpbypz3_eo text file with the staggered field for the ref. code. Returning...";   
        return;
	// */
	
	in.load(sf_in);
	in2.load(sf_in2);
	in3.load(sf_in3);
	alpha.load(&alpha_host);
	beta.load(&beta_host);

	logger.info() << "Run kernel";
	if(switcher)
	  device->saxpbypz_eoprec_device(&in, &in2, &in3, &alpha, &beta, &out);
	else
	  device->saxpbypz_eoprec_device(&in, &in2, &in3, alpha_host, beta_host, &out);

	logger.info() << "result:";
	hmc_float cpu_res;
	device->set_float_to_global_squarenorm_eoprec_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_sf_gaussian_staggered_eo(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "set_gaussian_spinorfield_stagg_eoprec";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);

	physics::PRNG prng(system);
	auto device = system.get_devices().at(0)->get_spinor_staggered_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_eoprec_spinorfieldsize(params);
	const SU3vec out(NUM_ELEMENTS_SF, device->get_device());

	int iterations = params.get_integrationsteps(0);

	su3vec * sf_out;
	sf_out = new su3vec[NUM_ELEMENTS_SF * iterations];
	BOOST_REQUIRE(sf_out);

	auto prng_buf = prng.get_buffers().at(0);

	hmc_float sum = 0;
	for (int i = 0; i< iterations; i++){
	  if(i%200==0)logger.info() << "Run kernel for the " << i << "th time";
	  device->set_gaussian_spinorfield_eoprec_device(&out, prng_buf);
	  out.dump(&sf_out[i*NUM_ELEMENTS_SF]);
	  //Here we sum the entries to calculate the mean later
	  sum += count_sf(&sf_out[i*NUM_ELEMENTS_SF], NUM_ELEMENTS_SF);
	}
	logger.info() << "result: mean";
	hmc_float cpu_res = 0.;
	//sum is the sum of iterations*NUM_ELEMENTS_SF*6 real gaussian numbers
	sum = sum/iterations/NUM_ELEMENTS_SF/6;	
	cpu_res= sum;
	logger.info() << cpu_res;

	if(params.get_read_multiple_configs()  == false){
	  //CP: calc std derivation
	  hmc_float var=0.;
	  for (int i=0; i<iterations; i++){
	    var += calc_var_sf(&sf_out[i*NUM_ELEMENTS_SF], NUM_ELEMENTS_SF, sum);
	  }
	  //var is the sum of iterations*NUM_ELEMENTS_SF*6 square deviations
	  var=var/iterations/NUM_ELEMENTS_SF/6;
	  
	  cpu_res = sqrt(var);
	  logger.info() << "result: variance";
	  logger.info() << cpu_res;
	}

	//The Kolmogorov_Smirnov test requires set of n samples with n around 1000
	//(to big n and to small n are not good choices for this test)
	//So we can consider each kernel result as a set (NUM_ELEMENTS_SF*6=1536 for a 4^4 lattice)
	vector<vector<hmc_float>> samples;
	vector<hmc_float> tmp;
	vector<hmc_float> tmp2;
	for(int i=0; i<iterations; i++){
	  vector<hmc_float> tmp;
	  for(uint j=0; j<NUM_ELEMENTS_SF; j++){
	    tmp2=reals_from_su3vec(sf_out[i*NUM_ELEMENTS_SF+j]);
	    tmp.insert(tmp.end(),tmp2.begin(),tmp2.end());
	    tmp2.clear();
	  }
	  samples.push_back(tmp);
	  tmp.clear();
	}
	logger.info() << "Running Kolmogorov_Smirnov test (it should take half a minute)...";
	logger.info() << "Kolmogorov_Smirnov frequency (of K+): " << std::setprecision(16) << Kolmogorov_Smirnov(samples,0.,sqrt(0.5)) << " ---> It should be close 0.98";
	
	if(params.get_read_multiple_configs()==true){
	  //Let us use the same sets of samples to make the mean test to 1,2 and 3 sigma
	  //Note that in the test BOOST_CHECK is used.
	  mean_test_multiple_set(samples,2.,0.,sqrt(0.5));
	  mean_test_multiple_set(samples,3.,0.,sqrt(0.5));
	  mean_test_multiple_set(samples,4.,0.,sqrt(0.5));
	}else{
	  //Let us use the same sets of samples to make the variance test to 1,2 and 3 sigma
	  //Note that in the test BOOST_CHECK is used.
	  variance_test_multiple_set(samples,2.,sqrt(0.5));
	  variance_test_multiple_set(samples,3.,sqrt(0.5));
	  variance_test_multiple_set(samples,4.,sqrt(0.5));
	}
	
	testFloatSizeAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");

}


void test_sf_sax_vectorized_and_squarenorm_staggered_eo(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "sax_vectorized_and_squarenorm_eoprec";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto device = system.get_devices().at(0)->get_spinor_staggered_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_eoprec_spinorfieldsize(params);
	int NUM_EQS = params.get_md_approx_ord();
	const SU3vec in(NUM_ELEMENTS_SF*NUM_EQS, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(NUM_EQS, device->get_device());
	hardware::buffers::Plain<hmc_float> alpha(NUM_EQS, device->get_device());

	std::vector<hmc_float> alpha_host(NUM_EQS, params.get_beta());
	logger.info() << "Using:";
	for(uint i=0; i<alpha_host.size(); i++){
	  alpha_host[i] += i*params.get_rho();
	  logger.info() << "  alpha[" << i << "] = " << alpha_host[i];
	}

	su3vec * sf_in;
	sf_in = new su3vec[NUM_ELEMENTS_SF * NUM_EQS];
	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) {
	  fill_sf_with_one(sf_in, NUM_ELEMENTS_SF * NUM_EQS);
	}
	else {
	  fill_sf_with_random(sf_in, NUM_ELEMENTS_SF * NUM_EQS, 123);
	}
	BOOST_REQUIRE(sf_in);

	//The following three lines are to be used to produce the ref_vec file needed to get the ref_value
        //---> Comment them out when the reference values have been obtained! 
        /*
        print_staggeredfield_eo_to_textfile("ref_vec_sax_eo",sf_in,params); 
        logger.info() << "Produced the ref_vec text file with the staggered field for the ref. code. Returning...";   
        return;
	// */
	
	in.load(sf_in);
	alpha.load(&alpha_host[0]);

	logger.info() << "Run kernel";
	device->sax_vectorized_and_squarenorm_eoprec_device(&in, &alpha, NUM_EQS, &sqnorm);
	
	logger.info() << "result:";
	std::vector<hmc_float> cpu_res(NUM_EQS);
	sqnorm.dump(&cpu_res[0]);
	for(uint i=0; i<cpu_res.size(); i++)
	  logger.info() << "  cpu_res[" << i << "] = " << cpu_res[i];

	testFloatAgainstInputparameters(std::accumulate(cpu_res.begin(), cpu_res.end(), 0.0), params);
	BOOST_MESSAGE("Test done");
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////


BOOST_AUTO_TEST_SUITE(BUILD)

BOOST_AUTO_TEST_CASE( BUILD_1 )
{
  test_build("/opencl_module_spinors_staggered_build_input_1");
}

BOOST_AUTO_TEST_CASE( BUILD_2 )
{
  test_build("/opencl_module_spinors_staggered_build_input_2");
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(SF_SQUARENORM)

BOOST_AUTO_TEST_CASE( SF_SQUARENORM_1 )
{
  test_sf_squarenorm_staggered("/sf_squarenorm_staggered_input_1");
}

BOOST_AUTO_TEST_CASE( SF_SQUARENORM_2 )
{
  test_sf_squarenorm_staggered("/sf_squarenorm_staggered_input_2");
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(SF_SQUARENORM_REDUCTION)

BOOST_AUTO_TEST_CASE( SF_SQUARENORM_REDUCTION_1 )
{
  test_sf_squarenorm_staggered("/sf_squarenorm_staggered_reduction_input_1");
}

BOOST_AUTO_TEST_CASE( SF_SQUARENORM_REDUCTION_2 )
{
  test_sf_squarenorm_staggered("/sf_squarenorm_staggered_reduction_input_2");
}

BOOST_AUTO_TEST_CASE( SF_SQUARENORM_REDUCTION_3 )
{
  test_sf_squarenorm_staggered("/sf_squarenorm_staggered_reduction_input_3");
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(SF_SCALAR_PRODUCT)

BOOST_AUTO_TEST_CASE( SF_SCALAR_PRODUCT_1 )
{
  test_sf_scalar_product_staggered("/sf_scalar_product_staggered_input_1");
}

BOOST_AUTO_TEST_CASE( SF_SCALAR_PRODUCT_2 )
{
  test_sf_scalar_product_staggered("/sf_scalar_product_staggered_input_2");
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(SF_SCALAR_PRODUCT_REDUCTION)

BOOST_AUTO_TEST_CASE( SF_SCALAR_PRODUCT_REDUCTION_1 )
{
  test_sf_scalar_product_staggered("/sf_scalar_product_staggered_reduction_input_1");
}

BOOST_AUTO_TEST_CASE( SF_SCALAR_PRODUCT_REDUCTION_2 )
{
  test_sf_scalar_product_staggered("/sf_scalar_product_staggered_reduction_input_2");
}

BOOST_AUTO_TEST_CASE( SF_SCALAR_PRODUCT_REDUCTION_3 )
{
  test_sf_scalar_product_staggered("/sf_scalar_product_staggered_reduction_input_3");
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(SF_ZERO)

BOOST_AUTO_TEST_CASE( SF_ZERO_1 )
{
  test_sf_cold("/sf_set_zero_staggered_input_1", false);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(SF_COLD)

BOOST_AUTO_TEST_CASE( SF_COLD_1 )
{
  test_sf_cold("/sf_set_cold_staggered_input_1", true);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(SF_SAX)

BOOST_AUTO_TEST_CASE( SF_SAX_1 )
{
  test_sf_sax_staggered("/sf_sax_staggered_input_1");
}

BOOST_AUTO_TEST_CASE( SF_SAX_2 )
{
  test_sf_sax_staggered("/sf_sax_staggered_input_2");
}

BOOST_AUTO_TEST_CASE( SF_SAX_3 )
{
  test_sf_sax_staggered("/sf_sax_staggered_input_3");
}

BOOST_AUTO_TEST_CASE( SF_SAX_4 )
{
  test_sf_sax_staggered("/sf_sax_staggered_input_4");
}

BOOST_AUTO_TEST_CASE( SF_SAX_5 )
{
  test_sf_sax_staggered("/sf_sax_staggered_input_5");
}

BOOST_AUTO_TEST_CASE( SF_SAX_6 )
{
  test_sf_sax_staggered("/sf_sax_staggered_input_6");
}

BOOST_AUTO_TEST_CASE( SF_SAX_7 )
{
  test_sf_sax_staggered("/sf_sax_staggered_input_7");
}

BOOST_AUTO_TEST_CASE( SF_SAX_8 )
{
  test_sf_sax_staggered("/sf_sax_staggered_input_8");
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(SF_SAXPY)

BOOST_AUTO_TEST_CASE( SF_SAXPY_1 )
{
  test_sf_saxpy_staggered("/sf_saxpy_staggered_input_1", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_2 )
{
  test_sf_saxpy_staggered("/sf_saxpy_staggered_input_2", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_3 )
{
  test_sf_saxpy_staggered("/sf_saxpy_staggered_input_3", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_4 )
{
  test_sf_saxpy_staggered("/sf_saxpy_staggered_input_4", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_5 )
{
  test_sf_saxpy_staggered("/sf_saxpy_staggered_input_5", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_6 )
{
  test_sf_saxpy_staggered("/sf_saxpy_staggered_input_6", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_7 )
{
  test_sf_saxpy_staggered("/sf_saxpy_staggered_input_7", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_8 )
{
  test_sf_saxpy_staggered("/sf_saxpy_staggered_input_8", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_9 )
{
  test_sf_saxpy_staggered("/sf_saxpy_staggered_input_9", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_10 )
{
  test_sf_saxpy_staggered("/sf_saxpy_staggered_input_10", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_11 )
{
  test_sf_saxpy_staggered("/sf_saxpy_staggered_input_11", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_12 )
{
  test_sf_saxpy_staggered("/sf_saxpy_staggered_input_12", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_13 )
{
  test_sf_saxpy_staggered("/sf_saxpy_staggered_input_13", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_14 )
{
  test_sf_saxpy_staggered("/sf_saxpy_staggered_input_14", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_15 )
{
  test_sf_saxpy_staggered("/sf_saxpy_staggered_input_15", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_16 )
{
  test_sf_saxpy_staggered("/sf_saxpy_staggered_input_16", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_17 )
{
  test_sf_saxpy_staggered("/sf_saxpy_staggered_input_17", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_18 )
{
  test_sf_saxpy_staggered("/sf_saxpy_staggered_input_18", true);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(SF_SAXPBYPZ)

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_1 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_1");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_2 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_2");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_3 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_3");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_4 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_4");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_5 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_5");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_6 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_6");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_7 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_7");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_8 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_8");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_9 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_9");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_10 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_10");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_11 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_11");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_12 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_12");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_13 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_13");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_14 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_14");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_15 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_15");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_16 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_16");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_17 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_17");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_18 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_18");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_19 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_19");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_20 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_20");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_21 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_21");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_22 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_22");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_23 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_23");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_24 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_24");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_25 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_25");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_26 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_26");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_27 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_27");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_28 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_28");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_29 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_29");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_30 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_30");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_31 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_31");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_32 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_32");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_33 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_33");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_34 )
{
  test_sf_saxpbypz_staggered("/sf_saxpbypz_staggered_input_34");
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(SF_GAUSSIAN)

BOOST_AUTO_TEST_CASE( SF_GAUSSIAN_1 )
{
  test_sf_gaussian_staggered("/sf_gaussian_staggered_input_1");
}

BOOST_AUTO_TEST_CASE( SF_GAUSSIAN_2 )
{
  test_sf_gaussian_staggered("/sf_gaussian_staggered_input_2");
}

BOOST_AUTO_TEST_CASE( SF_GAUSSIAN_3 )
{
  test_sf_gaussian_staggered("/sf_gaussian_staggered_input_3");
}

BOOST_AUTO_TEST_CASE( SF_GAUSSIAN_4 )
{
  test_sf_gaussian_staggered("/sf_gaussian_staggered_input_4");
}

BOOST_AUTO_TEST_SUITE_END()

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////           EVEN-ODD PRECONDITIONING             //////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

BOOST_AUTO_TEST_SUITE(SF_CONVERT_EO)

BOOST_AUTO_TEST_CASE( SF_CONVERT_EO_1 )
{
  test_sf_convert_to_eo_staggered("/sf_convert_staggered_eo_input_1");
}

BOOST_AUTO_TEST_CASE( SF_CONVERT_EO_2 )
{
  test_sf_convert_to_eo_staggered("/sf_convert_staggered_eo_input_2");
}

BOOST_AUTO_TEST_CASE( SF_CONVERT_EO_3 )
{
  test_sf_convert_from_eo_staggered("/sf_convert_staggered_eo_input_1");
}

BOOST_AUTO_TEST_CASE( SF_CONVERT_EO_4 )
{
  test_sf_convert_from_eo_staggered("/sf_convert_staggered_eo_input_2");
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(SF_SQUARENORM_EO)

BOOST_AUTO_TEST_CASE( SF_SQUARENORM_EO_1 )
{
  test_sf_squarenorm_staggered_eo("/sf_squarenorm_staggered_eo_input_1");
}

BOOST_AUTO_TEST_CASE( SF_SQUARENORM_EO_2 )
{
  test_sf_squarenorm_staggered_eo("/sf_squarenorm_staggered_eo_input_2");
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(SF_SQUARENORM_EO_REDUCTION)

BOOST_AUTO_TEST_CASE( SF_SQUARENORM_EO_REDUCTION_1 )
{
  test_sf_squarenorm_staggered_eo("/sf_squarenorm_staggered_eo_reduction_input_1");
}

BOOST_AUTO_TEST_CASE( SF_SQUARENORM_EO_REDUCTION_2 )
{
  test_sf_squarenorm_staggered_eo("/sf_squarenorm_staggered_eo_reduction_input_2");
}

BOOST_AUTO_TEST_CASE( SF_SQUARENORM_EO_REDUCTION_3 )
{
  test_sf_squarenorm_staggered_eo("/sf_squarenorm_staggered_eo_reduction_input_3");
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(SF_SCALAR_PRODUCT_EO)

BOOST_AUTO_TEST_CASE( SF_SCALAR_PRODUCT_EO_1 )
{
  test_sf_scalar_product_staggered_eo("/sf_scalar_product_staggered_eo_input_1");
}

BOOST_AUTO_TEST_CASE( SF_SCALAR_PRODUCT_EO_2 )
{
  test_sf_scalar_product_staggered_eo("/sf_scalar_product_staggered_eo_input_2");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SCALAR_PRODUCT_REAL_EO)

BOOST_AUTO_TEST_CASE( SF_SCALAR_PRODUCT_REAL_EO_1 )
{
  test_sf_scalar_product_staggered_eo("/sf_scalar_product_staggered_eo_input_1", true);
}

BOOST_AUTO_TEST_CASE( SF_SCALAR_PRODUCT_REAL_EO_2 )
{
  test_sf_scalar_product_staggered_eo("/sf_scalar_product_staggered_eo_reduction_input_1", true);
}

BOOST_AUTO_TEST_CASE( SF_SCALAR_PRODUCT_REAL_EO_3 )
{
  test_sf_scalar_product_staggered_eo("/sf_scalar_product_staggered_eo_reduction_input_2", true);
}

BOOST_AUTO_TEST_CASE( SF_SCALAR_PRODUCT_REAL_EO_4 )
{
  test_sf_scalar_product_staggered_eo("/sf_scalar_product_staggered_eo_reduction_input_3", true);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SCALAR_PRODUCT_EO_REDUCTION)

BOOST_AUTO_TEST_CASE( SF_SCALAR_PRODUCT_EO_REDUCTION_1 )
{
  test_sf_scalar_product_staggered_eo("/sf_scalar_product_staggered_eo_reduction_input_1");
}

BOOST_AUTO_TEST_CASE( SF_SCALAR_PRODUCT_EO_REDUCTION_2 )
{
  test_sf_scalar_product_staggered_eo("/sf_scalar_product_staggered_eo_reduction_input_2");
}

BOOST_AUTO_TEST_CASE( SF_SCALAR_PRODUCT_EO_REDUCTION_3 )
{
  test_sf_scalar_product_staggered_eo("/sf_scalar_product_staggered_eo_reduction_input_3");
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(SF_COLD_EO)

BOOST_AUTO_TEST_CASE( SF_COLD_EO_1 )
{
	test_sf_cold_staggered_eo("/sf_set_cold_staggered_eo_input_1", true);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(SF_ZERO_EO)

BOOST_AUTO_TEST_CASE( SF_ZERO_EO_1 )
{
  test_sf_cold_staggered_eo("/sf_set_zero_staggered_eo_input_1",  false);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(SF_SAX_CPLX_EO)

BOOST_AUTO_TEST_CASE( SF_SAX_CPLX_EO_1 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_1", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAX_CPLX_EO_2 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_2", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAX_CPLX_EO_3 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_3", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAX_CPLX_EO_4 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_4", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAX_CPLX_EO_5 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_5", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAX_CPLX_EO_6 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_6", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAX_CPLX_EO_7 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_7", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAX_CPLX_EO_8 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_8", 0);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SAX_CPLX_ARG_EO)

BOOST_AUTO_TEST_CASE( SF_SAX_CPLX_ARG_EO_1 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_1", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAX_CPLX_ARG_EO_2 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_2", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAX_CPLX_ARG_EO_3 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_3", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAX_CPLX_ARG_EO_4 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_4", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAX_CPLX_ARG_EO_5 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_5", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAX_CPLX_ARG_EO_6 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_6", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAX_CPLX_ARG_EO_7 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_7", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAX_CPLX_ARG_EO_8 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_8", 1);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SAX_REAL_EO)

BOOST_AUTO_TEST_CASE( SF_SAX_REAL_EO_1 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_1", 2);
}

BOOST_AUTO_TEST_CASE( SF_SAX_REAL_EO_2 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_2", 2);
}

BOOST_AUTO_TEST_CASE( SF_SAX_REAL_EO_3 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_5", 2);
}

BOOST_AUTO_TEST_CASE( SF_SAX_REAL_EO_4 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_6", 2);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SAX_REAL_ARG_EO)

BOOST_AUTO_TEST_CASE( SF_SAX_REAL_ARG_EO_1 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_1", 3);
}

BOOST_AUTO_TEST_CASE( SF_SAX_REAL_ARG_EO_2 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_2", 3);
}

BOOST_AUTO_TEST_CASE( SF_SAX_REAL_ARG_EO_3 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_5", 3);
}

BOOST_AUTO_TEST_CASE( SF_SAX_REAL_ARG_EO_4 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_6", 3);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SAX_REAL_VEC_EO)

BOOST_AUTO_TEST_CASE( SF_SAX_REAL_VEC_EO_1 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_1", 4);
}

BOOST_AUTO_TEST_CASE( SF_SAX_REAL_VEC_EO_2 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_2", 4);
}

BOOST_AUTO_TEST_CASE( SF_SAX_REAL_VEC_EO_3 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_5", 4);
}

BOOST_AUTO_TEST_CASE( SF_SAX_REAL_VEC_EO_4 )
{
  test_sf_sax_staggered_eo("/sf_sax_staggered_eo_input_6", 4);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(SF_SAXPY_CPLX_EO)

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_EO_1 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_1", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_EO_2 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_2", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_EO_3 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_3", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_EO_4 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_4", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_EO_5 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_5", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_EO_6 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_6", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_EO_7 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_7", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_EO_8 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_8", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_EO_9 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_9", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_EO_10 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_10", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_EO_11 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_11", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_EO_12 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_12", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_EO_13 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_13", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_EO_14 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_14", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_EO_15 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_15", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_EO_16 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_16", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_EO_17 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_17", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_EO_18 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_18", 0);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SAXPY_CPLX_ARG_EO)

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_ARG_EO_1 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_1", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_ARG_EO_2 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_2", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_ARG_EO_3 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_3", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_ARG_EO_4 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_4", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_ARG_EO_5 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_5", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_ARG_EO_6 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_6", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_ARG_EO_7 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_7", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_ARG_EO_8 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_8", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_ARG_EO_9 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_9", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_ARG_EO_10 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_10", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_ARG_EO_11 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_11", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_ARG_EO_12 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_12", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_ARG_EO_13 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_13", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_ARG_EO_14 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_14", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_ARG_EO_15 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_15", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_ARG_EO_16 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_16", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_ARG_EO_17 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_17", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_CPLX_ARG_EO_18 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_18", 1);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SAXPY_REAL_EO)

BOOST_AUTO_TEST_CASE( SF_SAXPY_REAL_EO_1 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_1", 2);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_REAL_EO_2 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_2", 2);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_REAL_EO_3 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_6", 2);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_REAL_EO_4 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_10", 2);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_REAL_EO_5 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_11", 2);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_REAL_EO_6 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_15", 2);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SAXPY_REAL_ARG_EO)

BOOST_AUTO_TEST_CASE( SF_SAXPY_REAL_ARG_EO_1 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_1", 3);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_REAL_ARG_EO_2 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_2", 3);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_REAL_ARG_EO_3 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_6", 3);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_REAL_ARG_EO_4 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_10", 3);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_REAL_ARG_EO_5 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_11", 3);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_REAL_ARG_EO_6 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_15", 3);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SAXPY_REAL_VEC_EO)

BOOST_AUTO_TEST_CASE( SF_SAXPY_REAL_VEC_EO_1 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_1", 4);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_REAL_VEC_EO_2 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_2", 4);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_REAL_VEC_EO_3 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_6", 4);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_REAL_VEC_EO_4 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_10", 4);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_REAL_VEC_EO_5 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_11", 4);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_REAL_VEC_EO_6 )
{
  test_sf_saxpy_staggered_eo("/sf_saxpy_staggered_eo_input_15", 4);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(SF_SAXPBY_EO)

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_1 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_1", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_2 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_2", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_3 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_3", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_4 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_4", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_5 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_5", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_6 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_6", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_7 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_7", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_8 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_8", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_9 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_9", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_10 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_10", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_11 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_11", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_12 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_12", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_13 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_13", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_14 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_14", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_15 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_15", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_16 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_16", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_17 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_17", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_18 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_18", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_19 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_19", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_20 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_20", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_21 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_21", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_22 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_22", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_23 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_23", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_24 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_24", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_25 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_25", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_26 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_26", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_27 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_27", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_28 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_28", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_29 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_29", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_30 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_30", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_31 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_31", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_32 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_32", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_33 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_33", 0);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_EO_34 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_34", 0);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SAXPBY_ARG_EO)

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_1 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_1", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_2 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_2", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_3 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_3", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_4 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_4", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_5 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_5", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_6 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_6", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_7 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_7", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_8 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_8", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_9 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_9", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_10 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_10", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_11 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_11", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_12 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_12", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_13 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_13", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_14 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_14", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_15 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_15", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_16 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_16", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_17 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_17", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_18 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_18", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_19 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_19", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_20 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_20", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_21 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_21", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_22 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_22", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_23 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_23", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_24 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_24", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_25 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_25", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_26 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_26", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_27 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_27", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_28 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_28", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_29 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_29", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_30 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_30", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_31 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_31", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_32 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_32", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_33 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_33", 1);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_ARG_EO_34 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_34", 1);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SAXPBY_REAL_EO)

BOOST_AUTO_TEST_CASE( SF_SAXPBY_REAL_EO_1 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_1", 2);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_REAL_EO_2 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_3", 2);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_REAL_EO_3 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_9", 2);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_REAL_EO_4 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_11", 2);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_REAL_EO_5 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_18", 2);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_REAL_EO_6 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_20", 2);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_REAL_EO_7 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_26", 2);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_REAL_EO_8 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_28", 2);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SAXPBY_REAL_ARG_EO)

BOOST_AUTO_TEST_CASE( SF_SAXPBY_REAL_ARG_EO_1 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_1", 3);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_REAL_ARG_EO_2 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_3", 3);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_REAL_ARG_EO_3 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_9", 3);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_REAL_ARG_EO_4 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_11", 3);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_REAL_ARG_EO_5 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_18", 3);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_REAL_ARG_EO_6 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_20", 3);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_REAL_ARG_EO_7 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_26", 3);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_REAL_ARG_EO_8 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_28", 3);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SAXPBY_REAL_VEC_EO)

BOOST_AUTO_TEST_CASE( SF_SAXPBY_REAL_VEC_EO_1 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_1", 4);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_REAL_VEC_EO_2 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_3", 4);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_REAL_VEC_EO_3 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_9", 4);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_REAL_VEC_EO_4 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_11", 4);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_REAL_VEC_EO_5 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_18", 4);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_REAL_VEC_EO_6 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_20", 4);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_REAL_VEC_EO_7 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_26", 4);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBY_REAL_VEC_EO_8 )
{
  test_sf_saxpby_staggered_eo("/sf_saxpby_staggered_eo_input_28", 4);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(SF_SAXPBYPZ_EO)

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_1 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_1");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_2 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_2");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_3 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_3");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_4 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_4");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_5 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_5");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_6 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_6");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_7 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_7");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_8 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_8");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_9 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_9");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_10 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_10");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_11 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_11");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_12 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_12");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_13 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_13");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_14 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_14");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_15 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_15");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_16 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_16");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_17 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_17");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_18 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_18");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_19 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_19");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_20 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_20");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_21 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_21");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_22 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_22");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_23 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_23");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_24 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_24");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_25 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_25");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_26 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_26");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_27 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_27");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_28 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_28");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_29 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_29");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_30 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_30");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_31 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_31");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_32 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_32");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_33 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_33");
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_EO_34 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_34");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SAXPBYPZ_ARG_EO)

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_1 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_1", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_2 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_2", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_3 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_3", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_4 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_4", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_5 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_5", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_6 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_6", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_7 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_7", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_8 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_8", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_9 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_9", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_10 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_10", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_11 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_11", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_12 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_12", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_13 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_13", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_14 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_14", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_15 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_15", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_16 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_16", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_17 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_17", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_18 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_18", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_19 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_19", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_20 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_20", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_21 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_21", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_22 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_22", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_23 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_23", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_24 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_24", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_25 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_25", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_26 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_26", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_27 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_27", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_28 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_28", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_29 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_29", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_30 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_30", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_31 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_31", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_32 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_32", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_33 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_33", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPBYPZ_ARG_EO_34 )
{
  test_sf_saxpbypz_staggered_eo("/sf_saxpbypz_staggered_eo_input_34", false);
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(SF_GAUSSIAN_EO)

BOOST_AUTO_TEST_CASE( SF_GAUSSIAN_EO_1 )
{
  test_sf_gaussian_staggered_eo("/sf_gaussian_staggered_eo_input_1");
}

BOOST_AUTO_TEST_CASE( SF_GAUSSIAN_EO_2 )
{
  test_sf_gaussian_staggered_eo("/sf_gaussian_staggered_eo_input_2");
}

BOOST_AUTO_TEST_CASE( SF_GAUSSIAN_EO_3 )
{
  test_sf_gaussian_staggered_eo("/sf_gaussian_staggered_eo_input_3");
}

BOOST_AUTO_TEST_CASE( SF_GAUSSIAN_EO_4 )
{
  test_sf_gaussian_staggered_eo("/sf_gaussian_staggered_eo_input_4");
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(SF_SAX_VEC_AND_SQNORM)

BOOST_AUTO_TEST_CASE( SF_SAX_VEC_AND_SQNORM_1 )
{
  test_sf_sax_vectorized_and_squarenorm_staggered_eo("/sf_sax_vec_sqnorm_staggered_eo_input_1");
}

BOOST_AUTO_TEST_CASE( SF_SAX_VEC_AND_SQNORM_2 )
{
  test_sf_sax_vectorized_and_squarenorm_staggered_eo("/sf_sax_vec_sqnorm_staggered_eo_input_2");
}

BOOST_AUTO_TEST_CASE( SF_SAX_VEC_AND_SQNORM_3 )
{
  test_sf_sax_vectorized_and_squarenorm_staggered_eo("/sf_sax_vec_sqnorm_staggered_eo_input_3");
}

BOOST_AUTO_TEST_CASE( SF_SAX_VEC_AND_SQNORM_4 )
{
  test_sf_sax_vectorized_and_squarenorm_staggered_eo("/sf_sax_vec_sqnorm_staggered_eo_input_4");
}

BOOST_AUTO_TEST_CASE( SF_SAX_VEC_AND_SQNORM_5 )
{
  test_sf_sax_vectorized_and_squarenorm_staggered_eo("/sf_sax_vec_sqnorm_staggered_eo_input_5");
}

BOOST_AUTO_TEST_CASE( SF_SAX_VEC_AND_SQNORM_6 )
{
  test_sf_sax_vectorized_and_squarenorm_staggered_eo("/sf_sax_vec_sqnorm_staggered_eo_input_6");
}

BOOST_AUTO_TEST_CASE( SF_SAX_VEC_AND_SQNORM_7 )
{
  test_sf_sax_vectorized_and_squarenorm_staggered_eo("/sf_sax_vec_sqnorm_staggered_eo_input_7");
}

BOOST_AUTO_TEST_CASE( SF_SAX_VEC_AND_SQNORM_8 )
{
  test_sf_sax_vectorized_and_squarenorm_staggered_eo("/sf_sax_vec_sqnorm_staggered_eo_input_8");
}

BOOST_AUTO_TEST_CASE( SF_SAX_VEC_AND_SQNORM_9 )
{
  test_sf_sax_vectorized_and_squarenorm_staggered_eo("/sf_sax_vec_sqnorm_staggered_eo_input_9");
}

BOOST_AUTO_TEST_CASE( SF_SAX_VEC_AND_SQNORM_10 )
{
  test_sf_sax_vectorized_and_squarenorm_staggered_eo("/sf_sax_vec_sqnorm_staggered_eo_input_10");
}

BOOST_AUTO_TEST_SUITE_END()

