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
	TestGaugefield cpu(&system);
	cl_int err = CL_SUCCESS;
	hardware::code::Spinors * device = cpu.get_device();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = meta::get_spinorfieldsize(params);
	const Plain<spinor> in(NUM_ELEMENTS_SF, device->get_device());
	const Plain<spinor> in2(NUM_ELEMENTS_SF, device->get_device());
	const Plain<spinor> out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
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
	if (switcher)
	  device->saxpy_device(&in, &in2, &alpha, &out);
	else
	  device->saxpy_device(&in, &in2, alpha_host, &out);

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

void test_sf_saxsbypz(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "saxsbypz";
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
	const Plain<spinor> in3(NUM_ELEMENTS_SF, device->get_device());
	const Plain<spinor> out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	hardware::buffers::Plain<hmc_complex> alpha(1, device->get_device());
	hardware::buffers::Plain<hmc_complex> beta(1, device->get_device());
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

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
	TestGaugefield cpu(&system);
	cl_int err = CL_SUCCESS;
	hardware::code::Spinors * device = cpu.get_device();

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = meta::get_eoprec_spinorfieldsize(params);
	const Spinor in(NUM_ELEMENTS_SF, device->get_device());
	const Spinor in2(NUM_ELEMENTS_SF, device->get_device());
	const Spinor out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
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
	if (switcher)
	  device->saxpy_eoprec_device(&in, &in2, &alpha, &out);
	else
	  device->saxpy_eoprec_device(&in, &in2, alpha_host, &out);

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

void test_sf_saxsbypz_eo(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "saxsbypz_eo";
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
	const Spinor in3(NUM_ELEMENTS_SF, device->get_device());
	const Spinor out(NUM_ELEMENTS_SF, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	hardware::buffers::Plain<hmc_complex> alpha(1, device->get_device());
	hardware::buffers::Plain<hmc_complex> beta(1, device->get_device());
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

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
	logger.info() << "Finalize device";
	cpu.finalize();

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_cplx(std::string inputfile, int switcher)
{
  //switcher chooses between product and ratio and convert
	using namespace hardware::buffers;

	std::string kernelName;
	if (switcher == 0)
	  kernelName = "product";
	else if(switcher == 1)
	  kernelName = "ratio";
	else if (switcher == 2)
	  kernelName = "convert";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	cl_int err = CL_SUCCESS;
	hardware::code::Spinors * device = cpu.get_device();

	logger.info() << "Fill buffers...";
	hardware::buffers::Plain<hmc_complex> sqnorm(1, device->get_device());
	hardware::buffers::Plain<hmc_complex> alpha(1, device->get_device());
	hardware::buffers::Plain<hmc_complex> beta(1, device->get_device());
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	hmc_complex alpha_host = {params.get_beta(), params.get_rho()};
	logger.info() << "Use alpha = (" << alpha_host.re << ","<< alpha_host.im <<")";
	hmc_complex beta_host = {params.get_kappa(), params.get_mu()};
	logger.info() << "Use beta = (" << beta_host.re << ","<< beta_host.im <<")";

	alpha.load(&alpha_host);
	beta.load(&beta_host);

	logger.info() << "Run kernel";
	if(switcher == 0)
	  device->set_complex_to_product_device(&alpha, &beta, &sqnorm);
	else if (switcher ==1)
	  device->set_complex_to_ratio_device(&alpha, &beta, &sqnorm);
	if(switcher == 2){
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

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CPLX_CONVERT)

BOOST_AUTO_TEST_CASE( CPLX_CONVERT_1 )
{
  test_cplx("/cplx_convert_input_1", 2);
}

BOOST_AUTO_TEST_CASE( CPLX_CONVERT_2 )
{
  test_cplx("/cplx_convert_input_2", 2);
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
