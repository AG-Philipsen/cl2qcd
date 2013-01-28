#include "../gaugefield_hybrid.h"
#include "../meta/util.hpp"
#include "../host_random.h"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE_SPINORS
#include <boost/test/unit_test.hpp>

//some functionality
#include "test_util.h"

class TestGaugefield : public Gaugefield_hybrid {

public:
	TestGaugefield(const hardware::System * system) : Gaugefield_hybrid(system), prng(*system) {
		auto inputfile = system->get_inputparameters();
		init(1, inputfile.get_use_gpu() ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU, prng);
		meta::print_info_inverter("test program", inputfile);
	};

	virtual void init_tasks();
	virtual void finalize_opencl();

  hardware::code::Correlator * get_device();
  physics::PRNG* get_prng();

private:
	physics::PRNG prng;
};

void TestGaugefield::init_tasks()
{
	opencl_modules = new hardware::code::Opencl_Module* [get_num_tasks()];
	opencl_modules[0] = get_device_for_task(0)->get_correlator_code();
}

void TestGaugefield::finalize_opencl()
{
	Gaugefield_hybrid::finalize_opencl();
}

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

hardware::code::Correlator* TestGaugefield::get_device()
{
  return static_cast<hardware::code::Correlator*>(opencl_modules[0]);
}
physics::PRNG* TestGaugefield::get_prng()
{
	return &prng;
}

void test_build(std::string inputfile)
{
	logger.info() << "build opencl_module_correlators";
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	logger.info() << "Finalize device";
	cpu.finalize();
	BOOST_MESSAGE("Test done");
}

void test_src_volume(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "create_volume_source";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);

	physics::PRNG * prng = cpu.get_prng();
	cl_int err = CL_SUCCESS;
	hardware::code::Correlator * device = cpu.get_device();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = meta::get_spinorfieldsize(params);
	const Plain<spinor> out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	//CP: run the kernel a couple of times times
	int iterations = params.get_integrationsteps(0);

	spinor * sf_out;
	sf_out = new spinor[NUM_ELEMENTS_SF * iterations];
	BOOST_REQUIRE(sf_out);

	auto spinor_code = device->get_device()->get_spinor_code();
	auto prng_buf = prng->get_buffers().at(0);

	hmc_float sum = 0;
	for (int i = 0; i< iterations; i++){
	  logger.info() << "Run kernel";
	  device->create_volume_source_device(&out, prng_buf);
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

	logger.info() << "Finalize device";
	cpu.finalize();
	
	if(params.get_sourcecontent() == meta::Inputparameters::one){
	  testFloatAgainstInputparameters(cpu_res, params);
	} else{
	  testFloatSizeAgainstInputparameters(cpu_res, params);
	}
	BOOST_MESSAGE("Test done");
}

void test_src_zslice(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "create_zslice_source";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);

	physics::PRNG * prng = cpu.get_prng();
	cl_int err = CL_SUCCESS;
	hardware::code::Correlator * device = cpu.get_device();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = meta::get_spinorfieldsize(params);
	//CP: this source does have a weight only on one slice
	size_t NUM_ELEMENTS_SRC = meta::get_volspace(params);
	const Plain<spinor> out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	//CP: run the kernel a couple of times times
	int iterations = params.get_integrationsteps(0);

	spinor * sf_out;
	sf_out = new spinor[NUM_ELEMENTS_SF * iterations];
	BOOST_REQUIRE(sf_out);

	auto spinor_code = device->get_device()->get_spinor_code();
	auto prng_buf = prng->get_buffers().at(0);

	hmc_float sum = 0;
	for (int i = 0; i< iterations; i++){
	  logger.info() << "Run kernel";
	  device->create_zslice_source_device(&out, prng_buf, params.get_source_z());
	  out.dump(&sf_out[i*NUM_ELEMENTS_SF]);
	  sum += count_sf(&sf_out[i*NUM_ELEMENTS_SF], NUM_ELEMENTS_SF);
	}
	logger.info() << "result: mean";
	hmc_float cpu_res = 0.;
	sum = sum/iterations/NUM_ELEMENTS_SRC/24;	
	cpu_res= sum;
	logger.info() << cpu_res;

	if(params.get_read_multiple_configs()  == false){
	  //CP: calc std derivation
	  hmc_float var=0.;
	  for (int i=0; i<iterations; i++){
	    var += calc_var_sf(&sf_out[i*NUM_ELEMENTS_SF], NUM_ELEMENTS_SF, sum);
	  }
	  var=var/iterations/NUM_ELEMENTS_SRC/24;
	  
	  cpu_res = sqrt(var);
	  logger.info() << "result: variance";
	  logger.info() << cpu_res;
	}

	logger.info() << "Finalize device";
	cpu.finalize();
	
	if(params.get_sourcecontent() == meta::Inputparameters::one){
	  testFloatAgainstInputparameters(cpu_res, params);
	} else{
	  testFloatSizeAgainstInputparameters(cpu_res, params);
	}
	BOOST_MESSAGE("Test done");
}

void test_src_tslice(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "create_timeslice_source";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);

	physics::PRNG * prng = cpu.get_prng();
	cl_int err = CL_SUCCESS;
	hardware::code::Correlator * device = cpu.get_device();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = meta::get_spinorfieldsize(params);
	//CP: this source does have a weight only on one slice
	size_t NUM_ELEMENTS_SRC = params.get_ntime() * params.get_nspace() * params.get_nspace();
	const Plain<spinor> out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	//CP: run the kernel a couple of times times
	int iterations = params.get_integrationsteps(0);

	spinor * sf_out;
	sf_out = new spinor[NUM_ELEMENTS_SF * iterations];
	BOOST_REQUIRE(sf_out);

	auto spinor_code = device->get_device()->get_spinor_code();
	auto prng_buf = prng->get_buffers().at(0);

	hmc_float sum = 0;
	for (int i = 0; i< iterations; i++){
	  logger.info() << "Run kernel";
	  device->create_timeslice_source_device(&out, prng_buf, params.get_source_z());
	  out.dump(&sf_out[i*NUM_ELEMENTS_SF]);
	  sum += count_sf(&sf_out[i*NUM_ELEMENTS_SF], NUM_ELEMENTS_SF);
	}
	logger.info() << "result: mean";
	hmc_float cpu_res = 0.;
	sum = sum/iterations/NUM_ELEMENTS_SRC/24;	
	cpu_res= sum;
	logger.info() << cpu_res;

	if(params.get_read_multiple_configs()  == false){
	  //CP: calc std derivation
	  hmc_float var=0.;
	  for (int i=0; i<iterations; i++){
	    var += calc_var_sf(&sf_out[i*NUM_ELEMENTS_SF], NUM_ELEMENTS_SF, sum);
	  }
	  var=var/iterations/NUM_ELEMENTS_SRC/24;
	  
	  cpu_res = sqrt(var);
	  logger.info() << "result: variance";
	  logger.info() << cpu_res;
	}

	logger.info() << "Finalize device";
	cpu.finalize();
	
	if(params.get_sourcecontent() == meta::Inputparameters::one){
	  testFloatAgainstInputparameters(cpu_res, params);
	} else{
	  testFloatSizeAgainstInputparameters(cpu_res, params);
	}
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


BOOST_AUTO_TEST_SUITE(SRC_VOLUME)

BOOST_AUTO_TEST_CASE( SRC_VOLUME_1 )
{
  test_src_volume("/src_volume_input_1");
}

BOOST_AUTO_TEST_CASE( SRC_VOLUME_2 )
{
  test_src_volume("/src_volume_input_2");
}

BOOST_AUTO_TEST_CASE( SRC_VOLUME_3 )
{
  test_src_volume("/src_volume_input_3");
}

BOOST_AUTO_TEST_CASE( SRC_VOLUME_4 )
{
  test_src_volume("/src_volume_input_4");
}

BOOST_AUTO_TEST_CASE( SRC_VOLUME_5 )
{
  test_src_volume("/src_volume_input_5");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SRC_ZSLICE)

BOOST_AUTO_TEST_CASE( SRC_ZSLICE_1 )
{
  test_src_zslice("/src_zslice_input_1");
}

BOOST_AUTO_TEST_CASE( SRC_ZSLICE_2 )
{
  test_src_zslice("/src_zslice_input_2");
}

BOOST_AUTO_TEST_CASE( SRC_ZSLICE_3 )
{
  test_src_zslice("/src_zslice_input_3");
}

BOOST_AUTO_TEST_CASE( SRC_ZSLICE_4 )
{
  test_src_zslice("/src_zslice_input_4");
}

BOOST_AUTO_TEST_CASE( SRC_ZSLICE_5 )
{
  test_src_zslice("/src_zslice_input_5");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SRC_TSLICE)

BOOST_AUTO_TEST_CASE( SRC_TSLICE_1 )
{
  test_src_tslice("/src_tslice_input_1");
}

BOOST_AUTO_TEST_CASE( SRC_TSLICE_2 )
{
  test_src_tslice("/src_tslice_input_2");
}

BOOST_AUTO_TEST_CASE( SRC_TSLICE_3 )
{
  test_src_tslice("/src_tslice_input_3");
}

BOOST_AUTO_TEST_CASE( SRC_TSLICE_4 )
{
  test_src_tslice("/src_tslice_input_4");
}

BOOST_AUTO_TEST_CASE( SRC_TSLICE_5 )
{
  test_src_tslice("/src_tslice_input_5");
}

BOOST_AUTO_TEST_SUITE_END()

