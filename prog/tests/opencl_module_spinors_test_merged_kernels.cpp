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
		meta::print_info_hmc("test program", inputfile);
	};

	virtual void init_tasks();
	virtual void finalize_opencl();

	hardware::code::Spinors * get_device();
  physics::PRNG* get_prng();

private:
	physics::PRNG prng;
};

void TestGaugefield::init_tasks()
{
	opencl_modules = new hardware::code::Opencl_Module* [get_num_tasks()];
	opencl_modules[0] = get_device_for_task(0)->get_spinor_code();
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

hardware::code::Spinors* TestGaugefield::get_device()
{
	return static_cast<hardware::code::Spinors*>(opencl_modules[0]);
}
physics::PRNG* TestGaugefield::get_prng()
{
	return &prng;
}

void test_build(std::string inputfile)
{
	logger.info() << "build opencl_module_spinors";
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	logger.info() << "Finalize device";
	cpu.finalize();
	BOOST_MESSAGE("Test done");
}

void test_sf_saxpy_AND_squarenorm_eo(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "saxpy_ANS_squarenorm_eo";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	cl_int err = CL_SUCCESS;
	hardware::code::Spinors * device = cpu.get_device();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = meta::get_eoprec_spinorfieldsize(params);
	const Spinor in(NUM_ELEMENTS_SF, device->get_device());
	const Spinor in2(NUM_ELEMENTS_SF, device->get_device());
	const Spinor out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_complex> sqnorm(1, device->get_device());
	hardware::buffers::Plain<hmc_complex> alpha(1, device->get_device());
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

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
	device->saxpy_AND_squarenorm_eo_device(&in, &in2, &alpha, &out, &sqnorm);

	logger.info() << "result:";
	hmc_float cpu_res;
	//spinor_code->set_float_to_global_squarenorm_eoprec_device(&out, &sqnorm);
	//sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;
	logger.info() << "Finalize device";
	cpu.finalize();

	testFloatAgainstInputparameters(cpu_res, params);
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


BOOST_AUTO_TEST_SUITE(SF_SAXPY_AND_SQUARENORM_EO)

BOOST_AUTO_TEST_CASE( SF_SAXPY_AND_SQUARENORM_EO_1 )
{
  test_sf_saxpy_AND_squarenorm_eo("/sf_saxpy_AND_squarenorm_eo_input_1");
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_AND_SQUARENORM_EO_2 )
{
  test_sf_saxpy_AND_squarenorm_eo("/sf_saxpy_AND_squarenorm_eo_input_2");
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_AND_SQUARENORM_EO_3 )
{
  test_sf_saxpy_AND_squarenorm_eo("/sf_saxpy_AND_squarenorm_eo_input_3");
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_AND_SQUARENORM_EO_4 )
{
  test_sf_saxpy_AND_squarenorm_eo("/sf_saxpy_AND_squarenorm_eo_input_4");
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_AND_SQUARENORM_EO_5 )
{
  test_sf_saxpy_AND_squarenorm_eo("/sf_saxpy_AND_squarenorm_eo_input_5");
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_AND_SQUARENORM_EO_6 )
{
  test_sf_saxpy_AND_squarenorm_eo("/sf_saxpy_AND_squarenorm_eo_input_6");
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_AND_SQUARENORM_EO_7 )
{
  test_sf_saxpy_AND_squarenorm_eo("/sf_saxpy_AND_squarenorm_eo_input_7");
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_AND_SQUARENORM_EO_8 )
{
  test_sf_saxpy_AND_squarenorm_eo("/sf_saxpy_AND_squarenorm_eo_input_8");
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_AND_SQUARENORM_EO_9 )
{
  test_sf_saxpy_AND_squarenorm_eo("/sf_saxpy_AND_squarenorm_eo_input_9");
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_AND_SQUARENORM_EO_10 )
{
  test_sf_saxpy_AND_squarenorm_eo("/sf_saxpy_AND_squarenorm_eo_input_10");
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_AND_SQUARENORM_EO_11 )
{
  test_sf_saxpy_AND_squarenorm_eo("/sf_saxpy_AND_squarenorm_eo_input_11");
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_AND_SQUARENORM_EO_12 )
{
  test_sf_saxpy_AND_squarenorm_eo("/sf_saxpy_AND_squarenorm_eo_input_12");
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_AND_SQUARENORM_EO_13 )
{
  test_sf_saxpy_AND_squarenorm_eo("/sf_saxpy_AND_squarenorm_eo_input_13");
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_AND_SQUARENORM_EO_14 )
{
  test_sf_saxpy_AND_squarenorm_eo("/sf_saxpy_AND_squarenorm_eo_input_14");
}

BOOST_AUTO_TEST_SUITE_END()
