#include "../meta/util.hpp"
#include "../host_random.h"
#include "../physics/prng.hpp"
#include "../hardware/device.hpp"
#include "../hardware/code/spinors.hpp"
#include "../hardware/code/complex.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE_SPINORS
#include <boost/test/unit_test.hpp>

//some functionality
#include "test_util.h"

void fill_sf_with_one(spinor * sf_in, int size)
{
	for(int i = 0; i < size; ++i) {
		sf_in[i].e0.e0 = hmc_complex_one;
		sf_in[i].e0.e1 = hmc_complex_one;
		sf_in[i].e0.e2 = hmc_complex_one;
		sf_in[i].e1.e0 = hmc_complex_one;
		sf_in[i].e1.e1 = hmc_complex_one;
		sf_in[i].e1.e2 = hmc_complex_one;
		sf_in[i].e2.e0 = hmc_complex_one;
		sf_in[i].e2.e1 = hmc_complex_one;
		sf_in[i].e2.e2 = hmc_complex_one;
		sf_in[i].e3.e0 = hmc_complex_one;
		sf_in[i].e3.e1 = hmc_complex_one;
		sf_in[i].e3.e2 = hmc_complex_one;
	}
	return;
}

void fill_sf_with_zero(spinor * sf_in, int size)
{
	for(int i = 0; i < size; ++i) {
		sf_in[i].e0.e0 = hmc_complex_zero;
		sf_in[i].e0.e1 = hmc_complex_zero;
		sf_in[i].e0.e2 = hmc_complex_zero;
		sf_in[i].e1.e0 = hmc_complex_zero;
		sf_in[i].e1.e1 = hmc_complex_zero;
		sf_in[i].e1.e2 = hmc_complex_zero;
		sf_in[i].e2.e0 = hmc_complex_zero;
		sf_in[i].e2.e1 = hmc_complex_zero;
		sf_in[i].e2.e2 = hmc_complex_zero;
		sf_in[i].e3.e0 = hmc_complex_zero;
		sf_in[i].e3.e1 = hmc_complex_zero;
		sf_in[i].e3.e2 = hmc_complex_zero;
	}
	return;
}

void fill_sf_with_one_eo(spinor * sf_in, int size, bool eo, meta::Inputparameters &params)
{
  int ns = params.get_nspace();
  int nt = params.get_ntime();
  int x,y,z,t;
  for (x = 0; x<ns;x++){
    for (y = 0; y<ns;y++){
      for (z = 0; z<ns;z++){
	for (t = 0; t<nt;t++){
	  int coord[4];
	  coord[0] = t;
	  coord[1] = x;
	  coord[2] = y;
	  coord[3] = z;
	  int nspace =  get_nspace(coord, params);
	  int global_pos = get_global_pos(nspace, t, params);
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
	  
	  sf_in[global_pos].e0.e0 = content;
	  sf_in[global_pos].e0.e1 = content;
	  sf_in[global_pos].e0.e2 = content;
	  sf_in[global_pos].e1.e0 = content;
	  sf_in[global_pos].e1.e1 = content;
	  sf_in[global_pos].e1.e2 = content;
	  sf_in[global_pos].e2.e0 = content;
	  sf_in[global_pos].e2.e1 = content;
	  sf_in[global_pos].e2.e2 = content;
	  sf_in[global_pos].e3.e0 = content;
	  sf_in[global_pos].e3.e1 = content;
	  sf_in[global_pos].e3.e2 = content;
	}}}}
  return;
}

hmc_float count_sf_eo(spinor * sf_in, int size, bool eo, meta::Inputparameters &params)
{
  int ns = params.get_nspace();
  int nt = params.get_ntime();
  int x,y,z,t;
  hmc_float sum = 0.;
  for (x = 0; x<ns;x++){
    for (y = 0; y<ns;y++){
      for (z = 0; z<ns;z++){
	for (t = 0; t<nt;t++){
	  int coord[4];
	  coord[0] = t;
	  coord[1] = x;
	  coord[2] = y;
	  coord[3] = z;
	  int nspace =  get_nspace(coord, params);
	  int global_pos = get_global_pos(nspace, t, params);
	  if (global_pos > size)
	    break;
	  if (
	      ( eo ==true && (x+y+z+t) %2 == 0) ||
	      ( eo ==false &&  (x+y+z+t) %2 == 1 )
	      )
	    {
	      int i = global_pos;
	      sum +=
		sf_in[i].e0.e0.re+ sf_in[i].e0.e0.im 
		+sf_in[i].e0.e1.re+ sf_in[i].e0.e1.im 
		+sf_in[i].e0.e2.re+ sf_in[i].e0.e2.im 
		+sf_in[i].e1.e0.re+ sf_in[i].e0.e0.im 
		+sf_in[i].e1.e1.re+ sf_in[i].e0.e1.im 
		+sf_in[i].e1.e2.re+ sf_in[i].e0.e2.im 
		+sf_in[i].e2.e0.re+ sf_in[i].e0.e0.im 
		+sf_in[i].e2.e1.re+ sf_in[i].e0.e1.im 
		+sf_in[i].e2.e2.re+ sf_in[i].e0.e1.im 
		+sf_in[i].e3.e0.re+ sf_in[i].e0.e0.im 
		+sf_in[i].e3.e1.re+ sf_in[i].e0.e1.im 
		+sf_in[i].e3.e2.re+ sf_in[i].e0.e1.im;
	    }
	  else{
	    continue;
	  }
	}}}}
  return sum;
}

hmc_float count_sf(spinor * sf_in, int size)
{
  hmc_float sum = 0.;
  for (int i = 0; i<size;i++){
    sum +=
       sf_in[i].e0.e0.re+ sf_in[i].e0.e0.im 
      +sf_in[i].e0.e1.re+ sf_in[i].e0.e1.im 
      +sf_in[i].e0.e2.re+ sf_in[i].e0.e2.im 
      +sf_in[i].e1.e0.re+ sf_in[i].e1.e0.im 
      +sf_in[i].e1.e1.re+ sf_in[i].e1.e1.im 
      +sf_in[i].e1.e2.re+ sf_in[i].e1.e2.im 
      +sf_in[i].e2.e0.re+ sf_in[i].e2.e0.im 
      +sf_in[i].e2.e1.re+ sf_in[i].e2.e1.im 
      +sf_in[i].e2.e2.re+ sf_in[i].e2.e2.im 
      +sf_in[i].e3.e0.re+ sf_in[i].e3.e0.im 
      +sf_in[i].e3.e1.re+ sf_in[i].e3.e1.im 
      +sf_in[i].e3.e2.re+ sf_in[i].e3.e2.im;
  }
  return sum;
}

hmc_float calc_var(hmc_float in, hmc_float mean){
  return (in - mean) * (in - mean);
}

hmc_float calc_var_sf(spinor * sf_in, int size, hmc_float sum){
  hmc_float var = 0.;
  for(int k = 0; k<size; k++){
    var +=
      calc_var( sf_in[k].e0.e0.re , sum) 
      + calc_var( sf_in[k].e0.e0.im , sum) 
      + calc_var( sf_in[k].e0.e1.re , sum)
      + calc_var( sf_in[k].e0.e1.im , sum) 
      + calc_var( sf_in[k].e0.e2.re , sum) 
      + calc_var( sf_in[k].e0.e2.im , sum) 
      + calc_var( sf_in[k].e1.e0.re , sum) 
      + calc_var( sf_in[k].e1.e0.im , sum) 
      + calc_var( sf_in[k].e1.e1.re , sum) 
      + calc_var( sf_in[k].e1.e1.im , sum) 
      + calc_var( sf_in[k].e1.e2.re , sum) 
      + calc_var( sf_in[k].e1.e2.im , sum) 
      + calc_var( sf_in[k].e2.e0.re , sum)
      + calc_var( sf_in[k].e2.e0.im , sum) 
      + calc_var( sf_in[k].e2.e1.re , sum)
      + calc_var( sf_in[k].e2.e1.im , sum) 
      + calc_var( sf_in[k].e2.e2.re , sum)
      + calc_var( sf_in[k].e2.e2.im , sum) 
      + calc_var( sf_in[k].e3.e0.re , sum)
      + calc_var( sf_in[k].e3.e0.im , sum) 
      + calc_var( sf_in[k].e3.e1.re , sum)
      + calc_var( sf_in[k].e3.e1.im , sum) 
      + calc_var( sf_in[k].e3.e2.re , sum)
      + calc_var( sf_in[k].e3.e2.im , sum);
  }
  return var;
}


void fill_sf_with_random(spinor * sf_in, int size, int seed)
{
	prng_init(seed);
	for(int i = 0; i < size; ++i) {
		sf_in[i].e0.e0.re = prng_double();
		sf_in[i].e0.e1.re = prng_double();
		sf_in[i].e0.e2.re = prng_double();
		sf_in[i].e1.e0.re = prng_double();
		sf_in[i].e1.e1.re = prng_double();
		sf_in[i].e1.e2.re = prng_double();
		sf_in[i].e2.e0.re = prng_double();
		sf_in[i].e2.e1.re = prng_double();
		sf_in[i].e2.e2.re = prng_double();
		sf_in[i].e3.e0.re = prng_double();
		sf_in[i].e3.e1.re = prng_double();
		sf_in[i].e3.e2.re = prng_double();

		sf_in[i].e0.e0.im = prng_double();
		sf_in[i].e0.e1.im = prng_double();
		sf_in[i].e0.e2.im = prng_double();
		sf_in[i].e1.e0.im = prng_double();
		sf_in[i].e1.e1.im = prng_double();
		sf_in[i].e1.e2.im = prng_double();
		sf_in[i].e2.e0.im = prng_double();
		sf_in[i].e2.e1.im = prng_double();
		sf_in[i].e2.e2.im = prng_double();
		sf_in[i].e3.e0.im = prng_double();
		sf_in[i].e3.e1.im = prng_double();
		sf_in[i].e3.e2.im = prng_double();
	}
	return;
}

void fill_sf_with_random(spinor * sf_in, int size)
{
	fill_sf_with_random(sf_in, size, 123456);
}

void test_build(std::string inputfile)
{
	logger.info() << "build opencl_module_spinors";
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	for(auto device: system.get_devices()) {
		device->get_spinor_code();
	}
	BOOST_MESSAGE("Test done");
}

void test_sf_squarenorm(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "global_squarenorm";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto * device = system.get_devices().at(0)->get_spinor_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_spinorfieldsize(params);
	const Plain<spinor> in(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	spinor * sf_in;
	sf_in = new spinor[NUM_ELEMENTS_SF];
	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	else fill_sf_with_random(sf_in, NUM_ELEMENTS_SF);
	BOOST_REQUIRE(sf_in);

	in.load(sf_in);

	auto spinor_code = device->get_device()->get_spinor_code();

	logger.info() << "Run kernel";
	logger.info() << "result:";
	hmc_float cpu_res;
	spinor_code->set_float_to_global_squarenorm_device(&in, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}


void test_sf_squarenorm_eo(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "global_squarenorm_eoprec";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto * device = system.get_devices().at(0)->get_spinor_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_eoprec_spinorfieldsize(params);
	const Spinor in(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	spinor * sf_in;
	sf_in = new spinor[NUM_ELEMENTS_SF];
	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	else fill_sf_with_random(sf_in, NUM_ELEMENTS_SF);
	BOOST_REQUIRE(sf_in);

	in.load(sf_in);

	auto spinor_code = device->get_device()->get_spinor_code();

	logger.info() << "Run kernel";
	logger.info() << "result:";
	hmc_float cpu_res;
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_sf_scalar_product(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "scalar_product";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto * device = system.get_devices().at(0)->get_spinor_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_spinorfieldsize(params);
	const Plain<spinor> in(NUM_ELEMENTS_SF, device->get_device());
	const Plain<spinor> in2(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_complex> sqnorm(1, device->get_device());

	spinor * sf_in;
	sf_in = new spinor[NUM_ELEMENTS_SF];
	spinor * sf_in2;
	sf_in2 = new spinor[NUM_ELEMENTS_SF];
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

	in.load(sf_in);
	in2.load(sf_in2);

	auto spinor_code = device->get_device()->get_spinor_code();

	logger.info() << "Run kernel";
	logger.info() << "result:";
	hmc_complex cpu_res_tmp;
	spinor_code->set_complex_to_scalar_product_device(&in, &in2, &sqnorm);
	sqnorm.dump(&cpu_res_tmp);
	hmc_float cpu_res = cpu_res_tmp.re + cpu_res_tmp.im;
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_sf_scalar_product_eo(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "scalar_product";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto * device = system.get_devices().at(0)->get_spinor_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_eoprec_spinorfieldsize(params);
	const Spinor in(NUM_ELEMENTS_SF, device->get_device());
	const Spinor in2(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_complex> sqnorm(1, device->get_device());

	spinor * sf_in;
	sf_in = new spinor[NUM_ELEMENTS_SF];
	spinor * sf_in2;
	sf_in2 = new spinor[NUM_ELEMENTS_SF];
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

	in.load(sf_in);
	in2.load(sf_in2);

	auto spinor_code = device->get_device()->get_spinor_code();

	logger.info() << "Run kernel";
	logger.info() << "result:";
	hmc_complex cpu_res_tmp;
	spinor_code->set_complex_to_scalar_product_eoprec_device(&in, &in2, &sqnorm);
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
	  kernelName = "set_spinorfield_cold";
	else
	  kernelName = "set_spinorfield_zero";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto * device = system.get_devices().at(0)->get_spinor_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_spinorfieldsize(params);
	const Plain<spinor> in(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	auto spinor_code = device->get_device()->get_spinor_code();

	logger.info() << "Run kernel";
        if(switcher)
	  device->set_spinorfield_cold_device(&in);
        else
          device->set_zero_spinorfield_device(&in);

	logger.info() << "result:";
	hmc_float cpu_res;
	spinor_code->set_float_to_global_squarenorm_device(&in, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_sf_cold_eo(std::string inputfile, bool switcher)
{
  //switcher decides if the sf is set to cold or zero
	using namespace hardware::buffers;

	std::string kernelName;
	if(switcher)
	  kernelName = "set_spinorfield_cold_eo";
	else
	  kernelName = "set_spinorfield_zero_eo";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto * device = system.get_devices().at(0)->get_spinor_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_eoprec_spinorfieldsize(params);
	const Spinor in(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	auto spinor_code = device->get_device()->get_spinor_code();

	logger.info() << "Run kernel";
	if(switcher)
	  device->set_eoprec_spinorfield_cold_device(&in);
	else
	  device->set_zero_spinorfield_eoprec_device(&in);
	logger.info() << "result:";
	hmc_float cpu_res;
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_sf_sax(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "sax";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto * device = system.get_devices().at(0)->get_spinor_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_spinorfieldsize(params);
	const Plain<spinor> in(NUM_ELEMENTS_SF, device->get_device());
	const Plain<spinor> out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	hardware::buffers::Plain<hmc_complex> alpha(1, device->get_device());

	hmc_complex alpha_host = {params.get_beta(), params.get_rho()};
	logger.info() << "Use alpha = (" << alpha_host.re << ","<< alpha_host.im <<")";

	spinor * sf_in;
	sf_in = new spinor[NUM_ELEMENTS_SF];
	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) {
	  fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	}
	else {
	  fill_sf_with_random(sf_in, NUM_ELEMENTS_SF, 123);
	}
	BOOST_REQUIRE(sf_in);

	in.load(sf_in);
	alpha.load(&alpha_host);

	auto spinor_code = device->get_device()->get_spinor_code();

	logger.info() << "Run kernel";
	device->sax_device(&in, &alpha, &out);

	logger.info() << "result:";
	hmc_float cpu_res;
	spinor_code->set_float_to_global_squarenorm_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_sf_saxpy(std::string inputfile, bool switcher)
{
  //switcher chooses between saxpy and saxpy_arg kernel, which have the same functionality
	using namespace hardware::buffers;

	std::string kernelName;
	if( switcher)
	  kernelName = "saxpy";
	else
	  kernelName = "saxpy_arg";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto * device = system.get_devices().at(0)->get_spinor_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_spinorfieldsize(params);
	const Plain<spinor> in(NUM_ELEMENTS_SF, device->get_device());
	const Plain<spinor> in2(NUM_ELEMENTS_SF, device->get_device());
	const Plain<spinor> out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	hardware::buffers::Plain<hmc_complex> alpha(1, device->get_device());

	hmc_complex alpha_host = {params.get_beta(), params.get_rho()};
	logger.info() << "Use alpha = (" << alpha_host.re << ","<< alpha_host.im <<")";

	spinor * sf_in;
	spinor * sf_in2;
	sf_in = new spinor[NUM_ELEMENTS_SF];
	sf_in2 = new spinor[NUM_ELEMENTS_SF];
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

	in.load(sf_in);
	in2.load(sf_in2);
	alpha.load(&alpha_host);

	auto spinor_code = device->get_device()->get_spinor_code();

	logger.info() << "Run kernel";
	if (switcher)
	  device->saxpy_device(&in, &in2, &alpha, &out);
	else
	  device->saxpy_device(&in, &in2, alpha_host, &out);

	logger.info() << "result:";
	hmc_float cpu_res;
	spinor_code->set_float_to_global_squarenorm_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_sf_saxsbypz(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "saxsbypz";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto * device = system.get_devices().at(0)->get_spinor_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_spinorfieldsize(params);
	const Plain<spinor> in(NUM_ELEMENTS_SF, device->get_device());
	const Plain<spinor> in2(NUM_ELEMENTS_SF, device->get_device());
	const Plain<spinor> in3(NUM_ELEMENTS_SF, device->get_device());
	const Plain<spinor> out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	hardware::buffers::Plain<hmc_complex> alpha(1, device->get_device());
	hardware::buffers::Plain<hmc_complex> beta(1, device->get_device());

	hmc_complex alpha_host = {params.get_beta(), params.get_rho()};
	hmc_complex beta_host = {params.get_kappa(), params.get_mu()};
	logger.info() << "Use alpha = (" << alpha_host.re << ","<< alpha_host.im <<")";
	logger.info() << "Use beta = (" << beta_host.re << ","<< beta_host.im <<")";

	spinor * sf_in;
	spinor * sf_in2;
	spinor * sf_in3;
	sf_in = new spinor[NUM_ELEMENTS_SF];
	sf_in2 = new spinor[NUM_ELEMENTS_SF];
	sf_in3 = new spinor[NUM_ELEMENTS_SF];
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

	in.load(sf_in);
	in2.load(sf_in2);
	in3.load(sf_in3);
	alpha.load(&alpha_host);
	beta.load(&beta_host);

	auto spinor_code = device->get_device()->get_spinor_code();

	logger.info() << "Run kernel";
	device->saxsbypz_device(&in, &in2, &in3, &alpha, &beta, &out);

	logger.info() << "result:";
	hmc_float cpu_res;
	spinor_code->set_float_to_global_squarenorm_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_sf_sax_eo(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "sax_eo";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto * device = system.get_devices().at(0)->get_spinor_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_eoprec_spinorfieldsize(params);
	const Spinor in(NUM_ELEMENTS_SF, device->get_device());
	const Spinor out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	hardware::buffers::Plain<hmc_complex> alpha(1, device->get_device());

	hmc_complex alpha_host = {params.get_beta(), params.get_rho()};
	logger.info() << "Use alpha = (" << alpha_host.re << ","<< alpha_host.im <<")";

	spinor * sf_in;
	sf_in = new spinor[NUM_ELEMENTS_SF];
	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) {
	  fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	}
	else {
	  fill_sf_with_random(sf_in, NUM_ELEMENTS_SF, 123);
	}
	BOOST_REQUIRE(sf_in);

	in.load(sf_in);
	alpha.load(&alpha_host);

	auto spinor_code = device->get_device()->get_spinor_code();

	logger.info() << "Run kernel";
	device->sax_eoprec_device(&in, &alpha, &out);

	logger.info() << "result:";
	hmc_float cpu_res;
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_sf_saxpy_eo(std::string inputfile, bool switcher)
{
  //switcher chooses between saxpy and saxpy_arg kernel, which have the same functionality
	using namespace hardware::buffers;

	std::string kernelName;
	if( switcher)
	  kernelName = "saxpy_eo";
	else
	  kernelName = "saxpy_eo_arg";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto * device = system.get_devices().at(0)->get_spinor_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_eoprec_spinorfieldsize(params);
	const Spinor in(NUM_ELEMENTS_SF, device->get_device());
	const Spinor in2(NUM_ELEMENTS_SF, device->get_device());
	const Spinor out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	hardware::buffers::Plain<hmc_complex> alpha(1, device->get_device());

	hmc_complex alpha_host = {params.get_beta(), params.get_rho()};
	logger.info() << "Use alpha = (" << alpha_host.re << ","<< alpha_host.im <<")";

	spinor * sf_in;
	spinor * sf_in2;
	sf_in = new spinor[NUM_ELEMENTS_SF];
	sf_in2 = new spinor[NUM_ELEMENTS_SF];
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

	in.load(sf_in);
	in2.load(sf_in2);
	alpha.load(&alpha_host);

	auto spinor_code = device->get_device()->get_spinor_code();

	logger.info() << "Run kernel";
	if (switcher)
	  device->saxpy_eoprec_device(&in, &in2, &alpha, &out);
	else
	  device->saxpy_eoprec_device(&in, &in2, alpha_host, &out);

	logger.info() << "result:";
	hmc_float cpu_res;
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_sf_saxsbypz_eo(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "saxsbypz_eo";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto * device = system.get_devices().at(0)->get_spinor_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_eoprec_spinorfieldsize(params);
	const Spinor in(NUM_ELEMENTS_SF, device->get_device());
	const Spinor in2(NUM_ELEMENTS_SF, device->get_device());
	const Spinor in3(NUM_ELEMENTS_SF, device->get_device());
	const Spinor out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	hardware::buffers::Plain<hmc_complex> alpha(1, device->get_device());
	hardware::buffers::Plain<hmc_complex> beta(1, device->get_device());

	hmc_complex alpha_host = {params.get_beta(), params.get_rho()};
	hmc_complex beta_host = {params.get_kappa(), params.get_mu()};
	logger.info() << "Use alpha = (" << alpha_host.re << ","<< alpha_host.im <<")";
	logger.info() << "Use beta = (" << beta_host.re << ","<< beta_host.im <<")";

	spinor * sf_in;
	spinor * sf_in2;
	spinor * sf_in3;
	sf_in = new spinor[NUM_ELEMENTS_SF];
	sf_in2 = new spinor[NUM_ELEMENTS_SF];
	sf_in3 = new spinor[NUM_ELEMENTS_SF];
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

	in.load(sf_in);
	in2.load(sf_in2);
	in3.load(sf_in3);
	alpha.load(&alpha_host);
	beta.load(&beta_host);

	auto spinor_code = device->get_device()->get_spinor_code();

	logger.info() << "Run kernel";
	device->saxsbypz_eoprec_device(&in, &in2, &in3, &alpha, &beta, &out);

	logger.info() << "result:";
	hmc_float cpu_res;
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
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
	meta::Inputparameters params = create_parameters(inputfile);
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

void test_sf_convert_to_eo(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "convert_to_eo";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto * device = system.get_devices().at(0)->get_spinor_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF_EO = hardware::code::get_eoprec_spinorfieldsize(params);
	size_t NUM_ELEMENTS_SF = hardware::code::get_spinorfieldsize(params);
	const Plain<spinor> in(NUM_ELEMENTS_SF, device->get_device());
	const Spinor in2(NUM_ELEMENTS_SF_EO, device->get_device());
	const Spinor in3(NUM_ELEMENTS_SF_EO, device->get_device());
	const Plain<spinor> out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	spinor * sf_in;
	sf_in = new spinor[NUM_ELEMENTS_SF];
	if(params.get_read_multiple_configs() )
	  fill_sf_with_one_eo(sf_in, NUM_ELEMENTS_SF, true, params);
	else
	  fill_sf_with_one_eo(sf_in, NUM_ELEMENTS_SF, false, params);
	BOOST_REQUIRE(sf_in);

	in.load(sf_in);

	auto spinor_code = device->get_device()->get_spinor_code();

	logger.info() << "Run kernel";
	device->convert_to_eoprec_device(&in2, &in3, &in);

	logger.info() << "result:";
	hmc_float cpu_res;
	if(params.get_read_multiple_configs() ){
	  spinor_code->set_float_to_global_squarenorm_eoprec_device(&in3, &sqnorm);
	  sqnorm.dump(&cpu_res);
	  logger.info() << cpu_res;
	  //CP: this must be zero since only the even sites should be filled!
	  BOOST_REQUIRE_CLOSE(cpu_res, 0., params.get_solver_prec());
	  spinor_code->set_float_to_global_squarenorm_eoprec_device(&in2, &sqnorm);
	  sqnorm.dump(&cpu_res);
	  logger.info() << cpu_res;
	}
	else{
	  spinor_code->set_float_to_global_squarenorm_eoprec_device(&in2, &sqnorm);
	  sqnorm.dump(&cpu_res);
	  logger.info() << cpu_res;
	  //CP: this must be zero since only the odd sites should be filled!
	  BOOST_REQUIRE_CLOSE(cpu_res, 0., params.get_solver_prec());
	  spinor_code->set_float_to_global_squarenorm_eoprec_device(&in3, &sqnorm);
	  sqnorm.dump(&cpu_res);
	  logger.info() << cpu_res;
	}

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_sf_convert_from_eo(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "convert_from_eo";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	auto * device = system.get_devices().at(0)->get_spinor_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF_EO = hardware::code::get_eoprec_spinorfieldsize(params);
	size_t NUM_ELEMENTS_SF = hardware::code::get_spinorfieldsize(params);
	const Spinor in2(NUM_ELEMENTS_SF_EO, device->get_device());
	const Spinor in3(NUM_ELEMENTS_SF_EO, device->get_device());
	const Plain<spinor> out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	spinor * sf_eo1;
	spinor * sf_eo2;
	spinor * sf_out;
	sf_eo1 = new spinor[NUM_ELEMENTS_SF_EO];
	sf_eo2 = new spinor[NUM_ELEMENTS_SF_EO];
	sf_out = new spinor[NUM_ELEMENTS_SF];
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

	auto spinor_code = device->get_device()->get_spinor_code();

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

void test_sf_gaussian(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "generate_gaussian_spinorfield";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);

	physics::PRNG prng(system);
	auto * device = system.get_devices().at(0)->get_spinor_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_spinorfieldsize(params);
	const Plain<spinor> out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	//CP: run the kernel a couple of times times
	int iterations = params.get_integrationsteps(0);

	spinor * sf_out;
	sf_out = new spinor[NUM_ELEMENTS_SF * iterations];
	BOOST_REQUIRE(sf_out);

	auto spinor_code = device->get_device()->get_spinor_code();
	auto prng_buf = prng.get_buffers().at(0);

	hmc_float sum = 0;
	for (int i = 0; i< iterations; i++){
	  logger.info() << "Run kernel";
	  device->generate_gaussian_spinorfield_device(&out, prng_buf);
	  out.dump(&sf_out[i*NUM_ELEMENTS_SF]);
	  sum += count_sf(&sf_out[i*NUM_ELEMENTS_SF], NUM_ELEMENTS_SF);
	}
	logger.info() << "result: mean";
	hmc_float cpu_res = 0.;
	sum = sum/iterations/NUM_ELEMENTS_SF/24;	
	cpu_res= sum;
	logger.info() << cpu_res;

	if(params.get_read_multiple_configs()  == false){
	  //CP: calc std derivation
	  hmc_float var=0.;
	  for (int i=0; i<iterations; i++){
	    var += calc_var_sf(&sf_out[i*NUM_ELEMENTS_SF], NUM_ELEMENTS_SF, sum);
	  }
	  var=var/iterations/NUM_ELEMENTS_SF/24;
	  
	  cpu_res = sqrt(var);
	  logger.info() << "result: variance";
	  logger.info() << cpu_res;
	}


	testFloatSizeAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");

}

void test_sf_gaussian_eo(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "generate_gaussian_spinorfield_eo";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);

	physics::PRNG prng(system);
	auto * device = system.get_devices().at(0)->get_spinor_code();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = hardware::code::get_eoprec_spinorfieldsize(params);
	const Spinor out(NUM_ELEMENTS_SF, device->get_device());

	//CP: run the kernel a couple of times times
	int iterations = params.get_integrationsteps(0);

	spinor * sf_out;
	sf_out = new spinor[NUM_ELEMENTS_SF * iterations];
	BOOST_REQUIRE(sf_out);

	auto spinor_code = device->get_device()->get_spinor_code();
	auto prng_buf = prng.get_buffers().at(0);

	hmc_float sum = 0;
	for (int i = 0; i< iterations; i++){
	  logger.info() << "Run kernel";
	  device->generate_gaussian_spinorfield_eo_device(&out, prng_buf);
	  out.dump(&sf_out[i*NUM_ELEMENTS_SF]);
	  sum += count_sf(&sf_out[i*NUM_ELEMENTS_SF], NUM_ELEMENTS_SF);
	}
	logger.info() << "result: mean";
	hmc_float cpu_res = 0.;
	sum = sum/iterations/NUM_ELEMENTS_SF/24;	
	cpu_res= sum;
	logger.info() << cpu_res;

	if(params.get_read_multiple_configs()  == false){
	  //CP: calc std derivation
	  hmc_float var=0.;
	  for (int i=0; i<iterations; i++){
	    var += calc_var_sf(&sf_out[i*NUM_ELEMENTS_SF], NUM_ELEMENTS_SF, sum);
	  }
	  var=var/iterations/NUM_ELEMENTS_SF/24;
	  
	  cpu_res = sqrt(var);
	  logger.info() << "result: variance";
	  logger.info() << cpu_res;
	}

	testFloatSizeAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");

}

BOOST_AUTO_TEST_SUITE(BUILD)

BOOST_AUTO_TEST_CASE( BUILD_1 )
{
  test_build("/opencl_module_spinors_build_input_1");
}

BOOST_AUTO_TEST_CASE( BUILD_2 )
{
	test_build("/opencl_module_spinors_build_input_2");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SQUARENORM_EO)

BOOST_AUTO_TEST_CASE( SF_SQUARENORM_EO_1 )
{
  test_sf_squarenorm_eo("/sf_squarenorm_eo_input_1");
}

BOOST_AUTO_TEST_CASE( SF_SQUARENORM_EO_2 )
{
  test_sf_squarenorm_eo("/sf_squarenorm_eo_input_2");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SQUARENORM_EO_REDUCTION)

BOOST_AUTO_TEST_CASE( SF_SQUARENORM_EO_REDUCTION_1 )
{
  test_sf_squarenorm_eo("/sf_squarenorm_eo_reduction_input_1");
  //test_sf_squarenorm_eo("/sf_squarenorm_eo_input_1");
}

BOOST_AUTO_TEST_CASE( SF_SQUARENORM_EO_REDUCTION_2 )
{
  test_sf_squarenorm_eo("/sf_squarenorm_eo_reduction_input_2");
}

BOOST_AUTO_TEST_CASE( SF_SQUARENORM_EO_REDUCTION_3 )
{
  test_sf_squarenorm_eo("/sf_squarenorm_eo_reduction_input_3");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SQUARENORM)

BOOST_AUTO_TEST_CASE( SF_SQUARENORM_1 )
{
  test_sf_squarenorm("/sf_squarenorm_input_1");
}

BOOST_AUTO_TEST_CASE( SF_SQUARENORM_2 )
{
  test_sf_squarenorm("/sf_squarenorm_input_2");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SQUARENORM_REDUCTION)

BOOST_AUTO_TEST_CASE( SF_SQUARENORM_REDUCTION_1 )
{
  test_sf_squarenorm("/sf_squarenorm_reduction_input_1");
}

BOOST_AUTO_TEST_CASE( SF_SQUARENORM_REDUCTION_2 )
{
  test_sf_squarenorm("/sf_squarenorm_reduction_input_2");
}

BOOST_AUTO_TEST_CASE( SF_SQUARENORM_REDUCTION_3 )
{
  test_sf_squarenorm("/sf_squarenorm_reduction_input_3");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SCALAR_PRODUCT)

BOOST_AUTO_TEST_CASE( SF_SCALAR_PRODUCT_1 )
{
  test_sf_scalar_product("/sf_scalar_product_input_1");
}

BOOST_AUTO_TEST_CASE( SF_SCALAR_PRODUCT_2 )
{
  test_sf_scalar_product("/sf_scalar_product_input_2");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SCALAR_PRODUCT_REDUCTION)

BOOST_AUTO_TEST_CASE( SF_SCALAR_PRODUCT_REDUCTION_1 )
{
  test_sf_scalar_product("/sf_scalar_product_reduction_input_1");
}

BOOST_AUTO_TEST_CASE( SF_SCALAR_PRODUCT_REDUCTION_2 )
{
  test_sf_scalar_product("/sf_scalar_product_reduction_input_2");
}

BOOST_AUTO_TEST_CASE( SF_SCALAR_PRODUCT_REDUCTION_3 )
{
  test_sf_scalar_product("/sf_scalar_product_reduction_input_3");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SCALAR_PRODUCT_EO)

BOOST_AUTO_TEST_CASE( SF_SCALAR_PRODUCT_EO_1 )
{
  test_sf_scalar_product_eo("/sf_scalar_product_eo_input_1");
}

BOOST_AUTO_TEST_CASE( SF_SCALAR_PRODUCT_EO_2 )
{
  test_sf_scalar_product_eo("/sf_scalar_product_eo_input_2");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SCALAR_PRODUCT_EO_REDUCTION)

BOOST_AUTO_TEST_CASE( SF_SCALAR_PRODUCT_EO_REDUCTION_1 )
{
  test_sf_scalar_product_eo("/sf_scalar_product_eo_reduction_input_1");
}

BOOST_AUTO_TEST_CASE( SF_SCALAR_PRODUCT_EO_REDUCTION_2 )
{
  test_sf_scalar_product_eo("/sf_scalar_product_eo_reduction_input_2");
}

BOOST_AUTO_TEST_CASE( SF_SCALAR_PRODUCT_EO_REDUCTION_3 )
{
  test_sf_scalar_product_eo("/sf_scalar_product_eo_reduction_input_3");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_COLD)

BOOST_AUTO_TEST_CASE( SF_COLD_1 )
{
  test_sf_cold("/sf_cold_input_1", true);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_COLD_EO)

BOOST_AUTO_TEST_CASE( SF_COLD_EO_1 )
{
  test_sf_cold_eo("/sf_cold_eo_input_1", true);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_ZERO)

BOOST_AUTO_TEST_CASE( SF_ZERO_1 )
{
  test_sf_cold("/sf_zero_input_1", false);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_ZERO_EO)

BOOST_AUTO_TEST_CASE( SF_ZERO_EO_1 )
{
  test_sf_cold_eo("/sf_zero_eo_input_1",  false);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SAX)

BOOST_AUTO_TEST_CASE( SF_SAX_1 )
{
  test_sf_sax("/sf_sax_input_1");
}

BOOST_AUTO_TEST_CASE( SF_SAX_2 )
{
  test_sf_sax("/sf_sax_input_2");
}

BOOST_AUTO_TEST_CASE( SF_SAX_3 )
{
  test_sf_sax("/sf_sax_input_3");
}

BOOST_AUTO_TEST_CASE( SF_SAX_4 )
{
  test_sf_sax("/sf_sax_input_4");
}

BOOST_AUTO_TEST_CASE( SF_SAX_5 )
{
  test_sf_sax("/sf_sax_input_5");
}

BOOST_AUTO_TEST_CASE( SF_SAX_6 )
{
  test_sf_sax("/sf_sax_input_6");
}

BOOST_AUTO_TEST_CASE( SF_SAX_7 )
{
  test_sf_sax("/sf_sax_input_7");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SAX_EO)

BOOST_AUTO_TEST_CASE( SF_SAX_EO_1 )
{
  test_sf_sax_eo("/sf_sax_eo_input_1");
}

BOOST_AUTO_TEST_CASE( SF_SAX_EO_2 )
{
  test_sf_sax_eo("/sf_sax_eo_input_2");
}

BOOST_AUTO_TEST_CASE( SF_SAX_EO_3 )
{
  test_sf_sax_eo("/sf_sax_eo_input_3");
}

BOOST_AUTO_TEST_CASE( SF_SAX_EO_4 )
{
  test_sf_sax_eo("/sf_sax_eo_input_4");
}

BOOST_AUTO_TEST_CASE( SF_SAX_EO_5 )
{
  test_sf_sax_eo("/sf_sax_eo_input_5");
}

BOOST_AUTO_TEST_CASE( SF_SAX_EO_6 )
{
  test_sf_sax_eo("/sf_sax_eo_input_6");
}

BOOST_AUTO_TEST_CASE( SF_SAX_EO_7 )
{
  test_sf_sax_eo("/sf_sax_eo_input_7");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SAXPY)

BOOST_AUTO_TEST_CASE( SF_SAXPY_1 )
{
  test_sf_saxpy("/sf_saxpy_input_1", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_2 )
{
  test_sf_saxpy("/sf_saxpy_input_2", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_3 )
{
  test_sf_saxpy("/sf_saxpy_input_3", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_4 )
{
  test_sf_saxpy("/sf_saxpy_input_4", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_5 )
{
  test_sf_saxpy("/sf_saxpy_input_5", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_6 )
{
  test_sf_saxpy("/sf_saxpy_input_6", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_7 )
{
  test_sf_saxpy("/sf_saxpy_input_7", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_8 )
{
  test_sf_saxpy("/sf_saxpy_input_8", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_9 )
{
  test_sf_saxpy("/sf_saxpy_input_9", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_10 )
{
  test_sf_saxpy("/sf_saxpy_input_10", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_11 )
{
  test_sf_saxpy("/sf_saxpy_input_11", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_12 )
{
  test_sf_saxpy("/sf_saxpy_input_12", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_13 )
{
  test_sf_saxpy("/sf_saxpy_input_13", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_14 )
{
  test_sf_saxpy("/sf_saxpy_input_14", true);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SAXPY_ARG)

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_1 )
{
  test_sf_saxpy("/sf_saxpy_input_1", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_2 )
{
  test_sf_saxpy("/sf_saxpy_input_2", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_3 )
{
  test_sf_saxpy("/sf_saxpy_input_3", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_4 )
{
  test_sf_saxpy("/sf_saxpy_input_4", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_5 )
{
  test_sf_saxpy("/sf_saxpy_input_5", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_6 )
{
  test_sf_saxpy("/sf_saxpy_input_6", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_7 )
{
  test_sf_saxpy("/sf_saxpy_input_7", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_8 )
{
  test_sf_saxpy("/sf_saxpy_input_8", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_9 )
{
  test_sf_saxpy("/sf_saxpy_input_9", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_10 )
{
  test_sf_saxpy("/sf_saxpy_input_10", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_11 )
{
  test_sf_saxpy("/sf_saxpy_input_11", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_12 )
{
  test_sf_saxpy("/sf_saxpy_input_12", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_13 )
{
  test_sf_saxpy("/sf_saxpy_input_13", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_14 )
{
  test_sf_saxpy("/sf_saxpy_input_14", false);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SAXPY_EO)

BOOST_AUTO_TEST_CASE( SF_SAXPY_EO_1 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_1", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_EO_2 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_2", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_EO_3 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_3", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_EO_4 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_4", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_EO_5 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_5", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_EO_6 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_6", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_EO_7 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_7", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_EO_8 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_8", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_EO_9 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_9", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_EO_10 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_10", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_EO_11 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_11", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_EO_12 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_12", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_EO_13 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_13", true);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_EO_14 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_14", true);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SAXPY_ARG_EO)

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_EO_1 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_1", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_EO_2 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_2", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_EO_3 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_3", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_EO_4 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_4", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_EO_5 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_5", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_EO_6 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_6", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_EO_7 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_7", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_EO_8 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_8", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_EO_9 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_9", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_EO_10 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_10", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_EO_11 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_11", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_EO_12 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_12", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_EO_13 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_13", false);
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_ARG_EO_14 )
{
  test_sf_saxpy_eo("/sf_saxpy_eo_input_14", false);
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

BOOST_AUTO_TEST_SUITE(SF_SAXSBYPZ)

BOOST_AUTO_TEST_CASE( SF_SAXSBYPZ_1 )
{
  test_sf_saxsbypz("/sf_saxsbypz_input_1");
}

BOOST_AUTO_TEST_CASE( SF_SAXSBYPZ_2 )
{
  test_sf_saxsbypz("/sf_saxsbypz_input_2");
}

BOOST_AUTO_TEST_CASE( SF_SAXSBYPZ_3 )
{
  test_sf_saxsbypz("/sf_saxsbypz_input_3");
}

BOOST_AUTO_TEST_CASE( SF_SAXSBYPZ_4 )
{
  test_sf_saxsbypz("/sf_saxsbypz_input_4");
}

BOOST_AUTO_TEST_CASE( SF_SAXSBYPZ_5 )
{
  test_sf_saxsbypz("/sf_saxsbypz_input_5");
}

BOOST_AUTO_TEST_CASE( SF_SAXSBYPZ_6 )
{
  test_sf_saxsbypz("/sf_saxsbypz_input_6");
}

BOOST_AUTO_TEST_CASE( SF_SAXSBYPZ_7 )
{
  test_sf_saxsbypz("/sf_saxsbypz_input_7");
}

BOOST_AUTO_TEST_CASE( SF_SAXSBYPZ_8 )
{
  test_sf_saxsbypz("/sf_saxsbypz_input_8");
}

BOOST_AUTO_TEST_CASE( SF_SAXSBYPZ_9 )
{
  test_sf_saxsbypz("/sf_saxsbypz_input_9");
}

BOOST_AUTO_TEST_CASE( SF_SAXSBYPZ_10 )
{
  test_sf_saxsbypz("/sf_saxsbypz_input_10");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_SAXSBYPZ_EO)

BOOST_AUTO_TEST_CASE( SF_SAXSBYPZ_EO_1 )
{
  test_sf_saxsbypz_eo("/sf_saxsbypz_eo_input_1");
}

BOOST_AUTO_TEST_CASE( SF_SAXSBYPZ_EO_2 )
{
  test_sf_saxsbypz_eo("/sf_saxsbypz_eo_input_2");
}

BOOST_AUTO_TEST_CASE( SF_SAXSBYPZ_EO_3 )
{
  test_sf_saxsbypz_eo("/sf_saxsbypz_eo_input_3");
}

BOOST_AUTO_TEST_CASE( SF_SAXSBYPZ_EO_4 )
{
  test_sf_saxsbypz_eo("/sf_saxsbypz_eo_input_4");
}

BOOST_AUTO_TEST_CASE( SF_SAXSBYPZ_EO_5 )
{
  test_sf_saxsbypz_eo("/sf_saxsbypz_eo_input_5");
}

BOOST_AUTO_TEST_CASE( SF_SAXSBYPZ_EO_6 )
{
  test_sf_saxsbypz_eo("/sf_saxsbypz_eo_input_6");
}

BOOST_AUTO_TEST_CASE( SF_SAXSBYPZ_EO_7 )
{
  test_sf_saxsbypz_eo("/sf_saxsbypz_eo_input_7");
}

BOOST_AUTO_TEST_CASE( SF_SAXSBYPZ_EO_8 )
{
  test_sf_saxsbypz_eo("/sf_saxsbypz_eo_input_8");
}

BOOST_AUTO_TEST_CASE( SF_SAXSBYPZ_EO_9 )
{
  test_sf_saxsbypz_eo("/sf_saxsbypz_eo_input_9");
}

BOOST_AUTO_TEST_CASE( SF_SAXSBYPZ_EO_10 )
{
  test_sf_saxsbypz_eo("/sf_saxsbypz_eo_input_10");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_CONVERT_EO)

BOOST_AUTO_TEST_CASE( SF_CONVERT_EO_1 )
{
  test_sf_convert_to_eo("/sf_convert_eo_input_1");
}

BOOST_AUTO_TEST_CASE( SF_CONVERT_EO_2 )
{
  test_sf_convert_to_eo("/sf_convert_eo_input_2");
}

BOOST_AUTO_TEST_CASE( SF_CONVERT_EO_3 )
{
  test_sf_convert_from_eo("/sf_convert_eo_input_1");
}

BOOST_AUTO_TEST_CASE( SF_CONVERT_EO_4 )
{
  test_sf_convert_from_eo("/sf_convert_eo_input_2");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_GAUSSIAN)

BOOST_AUTO_TEST_CASE( SF_GAUSSIAN_1 )
{
  test_sf_gaussian("/sf_gaussian_input_1");
}

BOOST_AUTO_TEST_CASE( SF_GAUSSIAN_2 )
{
  test_sf_gaussian("/sf_gaussian_input_2");
}

BOOST_AUTO_TEST_CASE( SF_GAUSSIAN_3 )
{
  test_sf_gaussian("/sf_gaussian_input_3");
}

BOOST_AUTO_TEST_CASE( SF_GAUSSIAN_4 )
{
  test_sf_gaussian("/sf_gaussian_input_4");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_GAUSSIAN_EO)

BOOST_AUTO_TEST_CASE( SF_GAUSSIAN_EO_1 )
{
  test_sf_gaussian_eo("/sf_gaussian_eo_input_1");
}

BOOST_AUTO_TEST_CASE( SF_GAUSSIAN_EO_2 )
{
  test_sf_gaussian_eo("/sf_gaussian_eo_input_2");
}

BOOST_AUTO_TEST_CASE( SF_GAUSSIAN_EO_3 )
{
  test_sf_gaussian_eo("/sf_gaussian_eo_input_3");
}

BOOST_AUTO_TEST_CASE( SF_GAUSSIAN_EO_4 )
{
  test_sf_gaussian_eo("/sf_gaussian_eo_input_4");
}

BOOST_AUTO_TEST_SUITE_END()
