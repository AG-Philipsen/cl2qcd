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

void test_sf_squarenorm(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "global_squarenorm";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	cl_int err = CL_SUCCESS;
	hardware::code::Spinors * device = cpu.get_device();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = meta::get_spinorfieldsize(params);
	const Plain<spinor> in(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	spinor * sf_in;
	sf_in = new spinor[NUM_ELEMENTS_SF];
	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	else fill_sf_with_random(sf_in, NUM_ELEMENTS_SF);
	BOOST_REQUIRE(sf_in);

	in.load(sf_in);

	auto spinor_code = device->get_device()->get_spinor_code();
	auto gf_code = device->get_device()->get_gaugefield_code();

	logger.info() << "Run kernel";
	logger.info() << "result:";
	hmc_float cpu_res;
	spinor_code->set_float_to_global_squarenorm_device(&in, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;
	logger.info() << "Finalize device";
	cpu.finalize();

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
	TestGaugefield cpu(&system);
	cl_int err = CL_SUCCESS;
	hardware::code::Spinors * device = cpu.get_device();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = meta::get_eoprec_spinorfieldsize(params);
	const Spinor in(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	spinor * sf_in;
	sf_in = new spinor[NUM_ELEMENTS_SF];
	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	else fill_sf_with_random(sf_in, NUM_ELEMENTS_SF);
	BOOST_REQUIRE(sf_in);

	in.load(sf_in);

	auto spinor_code = device->get_device()->get_spinor_code();
	auto gf_code = device->get_device()->get_gaugefield_code();

	logger.info() << "Run kernel";
	logger.info() << "result:";
	hmc_float cpu_res;
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;
	logger.info() << "Finalize device";
	cpu.finalize();

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
	TestGaugefield cpu(&system);
	cl_int err = CL_SUCCESS;
	hardware::code::Spinors * device = cpu.get_device();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = meta::get_spinorfieldsize(params);
	const Plain<spinor> in(NUM_ELEMENTS_SF, device->get_device());
	const Plain<spinor> in2(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_complex> sqnorm(1, device->get_device());
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

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
	auto gf_code = device->get_device()->get_gaugefield_code();

	logger.info() << "Run kernel";
	logger.info() << "result:";
	hmc_complex cpu_res_tmp;
	spinor_code->set_complex_to_scalar_product_device(&in, &in2, &sqnorm);
	sqnorm.dump(&cpu_res_tmp);
	hmc_float cpu_res = cpu_res_tmp.re + cpu_res_tmp.im;
	logger.info() << cpu_res;
	logger.info() << "Finalize device";
	cpu.finalize();

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
	TestGaugefield cpu(&system);
	cl_int err = CL_SUCCESS;
	hardware::code::Spinors * device = cpu.get_device();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = meta::get_eoprec_spinorfieldsize(params);
	const Spinor in(NUM_ELEMENTS_SF, device->get_device());
	const Spinor in2(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_complex> sqnorm(1, device->get_device());
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

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
	auto gf_code = device->get_device()->get_gaugefield_code();

	logger.info() << "Run kernel";
	logger.info() << "result:";
	hmc_complex cpu_res_tmp;
	spinor_code->set_complex_to_scalar_product_eoprec_device(&in, &in2, &sqnorm);
	sqnorm.dump(&cpu_res_tmp);
	hmc_float cpu_res = cpu_res_tmp.re + cpu_res_tmp.im;
	logger.info() << cpu_res;
	logger.info() << "Finalize device";
	cpu.finalize();

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
	TestGaugefield cpu(&system);
	cl_int err = CL_SUCCESS;
	hardware::code::Spinors * device = cpu.get_device();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = meta::get_spinorfieldsize(params);
	const Plain<spinor> in(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	auto spinor_code = device->get_device()->get_spinor_code();
	auto gf_code = device->get_device()->get_gaugefield_code();

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
	logger.info() << "Finalize device";
	cpu.finalize();

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
	TestGaugefield cpu(&system);
	cl_int err = CL_SUCCESS;
	hardware::code::Spinors * device = cpu.get_device();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = meta::get_eoprec_spinorfieldsize(params);
	const Spinor in(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	auto spinor_code = device->get_device()->get_spinor_code();
	auto gf_code = device->get_device()->get_gaugefield_code();

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
	logger.info() << "Finalize device";
	cpu.finalize();

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
	TestGaugefield cpu(&system);
	cl_int err = CL_SUCCESS;
	hardware::code::Spinors * device = cpu.get_device();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = meta::get_spinorfieldsize(params);
	const Plain<spinor> in(NUM_ELEMENTS_SF, device->get_device());
	const Plain<spinor> out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	hardware::buffers::Plain<hmc_complex> alpha(1, device->get_device());
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

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
	auto gf_code = device->get_device()->get_gaugefield_code();

	logger.info() << "Run kernel";
	device->sax_device(&in, &alpha, &out);

	logger.info() << "result:";
	hmc_float cpu_res;
	spinor_code->set_float_to_global_squarenorm_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;
	logger.info() << "Finalize device";
	cpu.finalize();

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
	TestGaugefield cpu(&system);
	cl_int err = CL_SUCCESS;
	hardware::code::Spinors * device = cpu.get_device();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = meta::get_eoprec_spinorfieldsize(params);
	const Spinor in(NUM_ELEMENTS_SF, device->get_device());
	const Spinor out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	hardware::buffers::Plain<hmc_complex> alpha(1, device->get_device());
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

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
	auto gf_code = device->get_device()->get_gaugefield_code();

	logger.info() << "Run kernel";
	device->sax_eoprec_device(&in, &alpha, &out);

	logger.info() << "result:";
	hmc_float cpu_res;
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
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


BOOST_AUTO_TEST_SUITE(SF_COLD)

BOOST_AUTO_TEST_CASE( SF_COLD_1 )
{
  test_sf_cold("/sf_cold_input_1", true);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SF_COLD_EO)

BOOST_AUTO_TEST_CASE( SF_COLD_EO_1 )
{
	test_sf_cold("/sf_cold_eo_input_1", true);
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
  test_sf_cold("/sf_zero_eo_input_1",  false);
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
