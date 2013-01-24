#include "../gaugefield_hybrid.h"
#include "../meta/util.hpp"
#include "../host_random.h"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE_FERMIONS
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

	hardware::code::Fermions * get_device();

private:
	physics::PRNG prng;
};

void TestGaugefield::init_tasks()
{
	opencl_modules = new hardware::code::Opencl_Module* [get_num_tasks()];
	opencl_modules[0] = get_device_for_task(0)->get_fermion_code();
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

hardware::code::Fermions* TestGaugefield::get_device()
{
	return static_cast<hardware::code::Fermions*>(opencl_modules[0]);
}

void test_build(std::string inputfile)
{
	logger.info() << "build opencl_module_hmc";
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	logger.info() << "Finalize device";
	cpu.finalize();
	BOOST_MESSAGE("Test done");
}

void test_m_fermion(std::string inputfile, int switcher)
{
	//switcher switches between similar functions
	//0: m_wilson (pure wilson)
	//1: m_tm_plus (twisted mass, upper flavour)
	//2: m_tm_minus (twisted mass, lower flavour)
	using namespace hardware::buffers;

	std::string kernelName;
	if(switcher == 0) {
		kernelName = "m_wilson";
	} else if(switcher == 1) {
		kernelName = "m_tm_plus";
	} else if(switcher == 2) {
		kernelName = "m_tm_minus";
	} else {
		logger.fatal() << "wrong parameter in test_m_fermion";
	}
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	cl_int err = CL_SUCCESS;
	hardware::code::Fermions * device = cpu.get_device();
	spinor * sf_in;
	spinor * sf_out;

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = meta::get_spinorfieldsize(params);

	sf_in = new spinor[NUM_ELEMENTS_SF];
	sf_out = new spinor[NUM_ELEMENTS_SF];

	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	else fill_sf_with_random(sf_in, NUM_ELEMENTS_SF);
	BOOST_REQUIRE(sf_in);

	const Plain<spinor> in(NUM_ELEMENTS_SF, device->get_device());
	in.load(sf_in);
	const Plain<spinor> out(NUM_ELEMENTS_SF, device->get_device());
	out.load(sf_in);
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	auto spinor_code = device->get_device()->get_spinor_code();
	auto gf_code = device->get_device()->get_gaugefield_code();

	logger.info() << "|phi|^2:";
	hmc_float cpu_back;
	spinor_code->set_float_to_global_squarenorm_device(&in, &sqnorm);
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;
	logger.info() << "Run kernel";
	if(switcher == 0) {
		device->M_wilson_device(&in, &out,  gf_code->get_gaugefield(), params.get_kappa());
	} else if(switcher == 1) {
		device->M_tm_plus_device(&in, &out,  gf_code->get_gaugefield(), params.get_kappa(), meta::get_mubar(params));
	} else if(switcher == 2) {
		device->M_tm_minus_device(&in, &out,  gf_code->get_gaugefield(), params.get_kappa(), meta::get_mubar(params));
	} else {
		logger.fatal() << "wrong parameter in test_m_fermion";
	}
	logger.info() << "result:";
	hmc_float cpu_res;
	spinor_code->set_float_to_global_squarenorm_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;
	logger.info() << "Finalize device";
	cpu.finalize();

	logger.info() << "Clear buffers";
	delete[] sf_in;
	delete[] sf_out;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_m_wilson(std::string inputfile)
{
	test_m_fermion(inputfile, 0);
}

void test_m_tm_plus(std::string inputfile)
{
	test_m_fermion(inputfile, 1);
}

void test_m_tm_minus(std::string inputfile)
{
	test_m_fermion(inputfile, 2);
}

hmc_float calc_sf_sum(size_t NUM_ELEMS, spinor * in)
{
	hmc_float res = 0.;
	for(int i = 0; i < NUM_ELEMS; i++) {
		spinor tmp = in[i];
		res +=
		  tmp.e0.e0.re + tmp.e0.e0.im +
		  tmp.e0.e1.re + tmp.e0.e1.im +
		  tmp.e0.e2.re + tmp.e0.e2.im +
		  tmp.e1.e0.re + tmp.e1.e0.im +
		  tmp.e1.e1.re + tmp.e1.e1.im +
		  tmp.e1.e2.re + tmp.e1.e2.im +
		  tmp.e2.e0.re + tmp.e2.e0.im +
		  tmp.e2.e1.re + tmp.e2.e1.im +
		  tmp.e3.e2.re + tmp.e2.e2.im +
		  tmp.e3.e0.re + tmp.e3.e0.im +
		  tmp.e3.e1.re + tmp.e3.e1.im +
		  tmp.e3.e2.re + tmp.e3.e2.im ;
	}
	return res;
}

void test_gamma5(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName = "gamma5";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	cl_int err = CL_SUCCESS;
	hardware::code::Fermions * device = cpu.get_device();
	spinor * sf_in;

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = meta::get_spinorfieldsize(params);

	sf_in = new spinor[NUM_ELEMENTS_SF];

	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	else fill_sf_with_random(sf_in, NUM_ELEMENTS_SF);
	BOOST_REQUIRE(sf_in);

	const Plain<spinor> in(NUM_ELEMENTS_SF, device->get_device());
	in.load(sf_in);
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	auto spinor_code = device->get_device()->get_spinor_code();

	logger.info() << "|phi|^2:";
	hmc_float cpu_back;
	spinor_code->set_float_to_global_squarenorm_device(&in, &sqnorm);
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;
	logger.info() << "Run kernel";
	device->gamma5_device(&in);
	in.dump(sf_in);
	logger.info() << "result:";
	hmc_float cpu_res;
	cpu_res = calc_sf_sum(NUM_ELEMENTS_SF, sf_in);
	logger.info() << cpu_res;
	logger.info() << "Finalize device";
	cpu.finalize();

	logger.info() << "Clear buffers";
	delete[] sf_in;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_gamma5_eo(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName = "gamma5_eo";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	cl_int err = CL_SUCCESS;
	hardware::code::Fermions * device = cpu.get_device();
	spinor * sf_in;

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = meta::get_eoprec_spinorfieldsize(params);

	sf_in = new spinor[NUM_ELEMENTS_SF];

	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	else fill_sf_with_random(sf_in, NUM_ELEMENTS_SF);
	BOOST_REQUIRE(sf_in);

	const Spinor in(NUM_ELEMENTS_SF, device->get_device());
	in.load(sf_in);
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	auto spinor_code = device->get_device()->get_spinor_code();

	logger.info() << "|phi|^2:";
	hmc_float cpu_back;
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in, &sqnorm);

	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;
	logger.info() << "Run kernel";
	device->gamma5_eo_device(&in);
	in.dump(sf_in);
	logger.info() << "result:";
	hmc_float cpu_res;
	cpu_res = calc_sf_sum(NUM_ELEMENTS_SF, sf_in);
	logger.info() << cpu_res;
	logger.info() << "Finalize device";
	cpu.finalize();

	logger.info() << "Clear buffers";
	delete[] sf_in;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_m_tm_sitediagonal_plus_minus(std::string inputfile, bool switcher)
{
	using namespace hardware::buffers;

	std::string kernelName;
	if(switcher)
		kernelName = "m_tm_sitediagonal";
	else
		kernelName = "m_tm_sitediagonal_minus";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	hardware::code::Fermions * device = cpu.get_device();
	spinor * sf_in;
	spinor * sf_out;

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = meta::get_eoprec_spinorfieldsize(params);

	sf_in = new spinor[NUM_ELEMENTS_SF];
	sf_out = new spinor[NUM_ELEMENTS_SF];

	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	else fill_sf_with_random(sf_in, NUM_ELEMENTS_SF);
	BOOST_REQUIRE(sf_in);

	const Spinor in(NUM_ELEMENTS_SF, device->get_device());
	const Spinor out(NUM_ELEMENTS_SF, device->get_device());
	in.load(sf_in);
	out.load(sf_in);
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	auto spinor_code = device->get_device()->get_spinor_code();

	logger.info() << "|phi|^2:";
	hmc_float cpu_back;
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in, &sqnorm);
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;

	//switch according to "use_pointsource"
	hmc_float cpu_res;
	if(switcher) {
		device->M_tm_sitediagonal_device( &in, &out);
	} else {
		device->M_tm_sitediagonal_minus_device( &in, &out);
	}

	spinor_code->set_float_to_global_squarenorm_eoprec_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << "result:";
	logger.info() << cpu_res;
	logger.info() << "Finalize device";
	cpu.finalize();

	logger.info() << "Clear buffers";
	delete[] sf_in;
	delete[] sf_out;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_m_tm_sitediagonal(std::string inputfile)
{
	test_m_tm_sitediagonal_plus_minus(inputfile, true);
}

void test_m_tm_sitediagonal_minus(std::string inputfile)
{
	test_m_tm_sitediagonal_plus_minus(inputfile, false);
}

void test_m_tm_inverse_sitediagonal_plus_minus(std::string inputfile, bool switcher)
{
	using namespace hardware::buffers;

	std::string kernelName;
	if(switcher)
		kernelName = "m_tm_inverse_sitediagonal";
	else
		kernelName = "m_tm_inverse_sitediagonal_minus";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	hardware::code::Fermions * device = cpu.get_device();
	spinor * sf_in;
	spinor * sf_out;

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = meta::get_eoprec_spinorfieldsize(params);

	sf_in = new spinor[NUM_ELEMENTS_SF];
	sf_out = new spinor[NUM_ELEMENTS_SF];

	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	else fill_sf_with_random(sf_in, NUM_ELEMENTS_SF);
	BOOST_REQUIRE(sf_in);

	const Spinor in(NUM_ELEMENTS_SF, device->get_device());
	const Spinor out(NUM_ELEMENTS_SF, device->get_device());
	in.load(sf_in);
	out.load(sf_in);
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	auto spinor_code = device->get_device()->get_spinor_code();

	logger.info() << "|phi|^2:";
	hmc_float cpu_back;
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in, &sqnorm);
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;

	//switch according to "use_pointsource"
	hmc_float cpu_res;
	if(switcher) {
		device->M_tm_inverse_sitediagonal_device( &in, &out);
	} else {
		device->M_tm_inverse_sitediagonal_minus_device( &in, &out);
	}
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << "result:";
	logger.info() << cpu_res;
	logger.info() << "Finalize device";
	cpu.finalize();

	logger.info() << "Clear buffers";
	delete[] sf_in;
	delete[] sf_out;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_m_tm_inverse_sitediagonal(std::string inputfile)
{
	test_m_tm_inverse_sitediagonal_plus_minus(inputfile, true);
}

void test_m_tm_inverse_sitediagonal_minus(std::string inputfile)
{
	test_m_tm_inverse_sitediagonal_plus_minus(inputfile, false);
}

void test_dslash_eo(std::string inputfile)
{
	using namespace hardware::buffers;
	std::string kernelName = "dslash_eo";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	hardware::code::Fermions * device = cpu.get_device();

	logger.info() << "Fill buffers...";
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	size_t NUM_ELEMENTS_SF_EO = meta::get_eoprec_spinorfieldsize(params);
	spinor * sf_in_eo;
	sf_in_eo = new spinor[NUM_ELEMENTS_SF_EO];
	const Spinor in_eo_even(NUM_ELEMENTS_SF_EO, device->get_device());
	const Spinor out_eo(NUM_ELEMENTS_SF_EO, device->get_device());
	if(params.get_solver() == meta::Inputparameters::cg) fill_sf_with_one(sf_in_eo, NUM_ELEMENTS_SF_EO);
	else fill_sf_with_random(sf_in_eo, NUM_ELEMENTS_SF_EO);
	in_eo_even.load(sf_in_eo);

	auto spinor_code = device->get_device()->get_spinor_code();
	auto gf_code = device->get_device()->get_gaugefield_code();

	logger.info() << "|phi|^2:";
	hmc_float cpu_back;
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in_eo_even, &sqnorm);
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;

	//switch according to "use_pointsource"
	hmc_float cpu_res;
	if(params.get_use_pointsource()) {
		device->dslash_eo_device( &in_eo_even, &out_eo, gf_code->get_gaugefield(), EVEN, params.get_kappa() );
	} else {
		device->dslash_eo_device( &in_eo_even, &out_eo, gf_code->get_gaugefield(), ODD, params.get_kappa() );
	}
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&out_eo, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << "result:";
	logger.info() << cpu_res;
	logger.info() << "Finalize device";
	cpu.finalize();

	logger.info() << "Clear buffers";
	delete[] sf_in_eo;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

BOOST_AUTO_TEST_SUITE(BUILD)

BOOST_AUTO_TEST_CASE( BUILD_1 )
{
	test_build("/opencl_module_fermions_build_input_1");
}

BOOST_AUTO_TEST_CASE( BUILD_2 )
{
	test_build("/opencl_module_fermions_build_input_2");
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( M_WILSON )

BOOST_AUTO_TEST_CASE( M_WILSON_1)
{
	test_m_wilson("/m_wilson_input_1");
}

BOOST_AUTO_TEST_CASE( M_WILSON_2)
{
	test_m_wilson("/m_wilson_input_2");
}

BOOST_AUTO_TEST_CASE( M_WILSON_3)
{
	test_m_wilson("/m_wilson_input_3");
}

BOOST_AUTO_TEST_CASE( M_WILSON_4)
{
	test_m_wilson("/m_wilson_input_4");
}

BOOST_AUTO_TEST_CASE( M_WILSON_5)
{
	test_m_wilson("/m_wilson_input_5");
}

BOOST_AUTO_TEST_CASE( M_WILSON_6)
{
	test_m_wilson("/m_wilson_input_6");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( M_TM_MINUS  )

BOOST_AUTO_TEST_CASE( M_TM_MINUS_1 )
{
	test_m_tm_minus("/m_tm_minus_input_1");
}

BOOST_AUTO_TEST_CASE( M_TM_MINUS_2 )
{
	test_m_tm_minus("/m_tm_minus_input_2");
}

BOOST_AUTO_TEST_CASE( M_TM_MINUS_3 )
{
	test_m_tm_minus("/m_tm_minus_input_3");
}

BOOST_AUTO_TEST_CASE( M_TM_MINUS_4 )
{
	test_m_tm_minus("/m_tm_minus_input_4");
}

BOOST_AUTO_TEST_CASE( M_TM_MINUS_5 )
{
	test_m_tm_minus("/m_tm_minus_input_5");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( M_TM_PLUS )

BOOST_AUTO_TEST_CASE( M_TM_PLUS_1 )
{
	test_m_tm_plus("/m_tm_plus_input_1");
}

BOOST_AUTO_TEST_CASE( M_TM_PLUS_2 )
{
	test_m_tm_plus("/m_tm_plus_input_2");
}

BOOST_AUTO_TEST_CASE( M_TM_PLUS_3 )
{
	test_m_tm_plus("/m_tm_plus_input_3");
}

BOOST_AUTO_TEST_CASE( M_TM_PLUS_4 )
{
	test_m_tm_plus("/m_tm_plus_input_4");
}

BOOST_AUTO_TEST_CASE( M_TM_PLUS_5 )
{
	test_m_tm_plus("/m_tm_plus_input_5");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( GAMMA5 )

BOOST_AUTO_TEST_CASE( GAMMA5_1)
{
	test_gamma5("/gamma5_input_1");
}

BOOST_AUTO_TEST_CASE( GAMMA5_2 )
{
	test_gamma5("/gamma5_input_2");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( GAMMA5_EO)

BOOST_AUTO_TEST_CASE( GAMMA5_EO_1)
{
	test_gamma5_eo("/gamma5_eo_input_1");
}

BOOST_AUTO_TEST_CASE( GAMMA5_EO_2 )
{
	test_gamma5_eo("/gamma5_eo_input_2");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(M_TM_SITEDIAGONAL )

BOOST_AUTO_TEST_CASE( M_TM_SITEDIAGONAL_1)
{
	test_m_tm_sitediagonal("/m_tm_sitediagonal_input_1");
}

BOOST_AUTO_TEST_CASE( M_TM_SITEDIAGONAL_2)
{
	test_m_tm_sitediagonal("/m_tm_sitediagonal_input_2");
}

BOOST_AUTO_TEST_CASE( M_TM_SITEDIAGONAL_3)
{
	test_m_tm_sitediagonal("/m_tm_sitediagonal_input_3");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( M_TM_INVERSE_SITEDIAGONAL )

BOOST_AUTO_TEST_CASE( M_TM_INVERSE_SITEDIAGONAL_1)
{
	test_m_tm_inverse_sitediagonal("/m_tm_inverse_sitediagonal_input_1");
}

BOOST_AUTO_TEST_CASE( M_TM_INVERSE_SITEDIAGONAL_2)
{
	test_m_tm_inverse_sitediagonal("/m_tm_inverse_sitediagonal_input_2");
}

BOOST_AUTO_TEST_CASE( M_TM_INVERSE_SITEDIAGONAL_3)
{
	test_m_tm_inverse_sitediagonal("/m_tm_inverse_sitediagonal_input_3");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(M_TM_SITEDIAGONAL_MINUS )

BOOST_AUTO_TEST_CASE( M_TM_SITEDIAGONAL_MINUS_1)
{
	test_m_tm_sitediagonal_minus("/m_tm_sitediagonal_minus_input_1");
}

BOOST_AUTO_TEST_CASE( M_TM_SITEDIAGONAL_MINUS_2)
{
	test_m_tm_sitediagonal_minus("/m_tm_sitediagonal_minus_input_2");
}

BOOST_AUTO_TEST_CASE( M_TM_SITEDIAGONAL_MINUS_3)
{
	test_m_tm_sitediagonal_minus("/m_tm_sitediagonal_minus_input_3");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( M_TM_INVERSE_SITEDIAGONAL_MINUS)

BOOST_AUTO_TEST_CASE( M_TM_INVERSE_SITEDIAGONAL_MINUS_1)
{
	test_m_tm_inverse_sitediagonal_minus("/m_tm_inverse_sitediagonal_minus_input_1");
}

BOOST_AUTO_TEST_CASE( M_TM_INVERSE_SITEDIAGONAL_MINUS_2)
{
	test_m_tm_inverse_sitediagonal_minus("/m_tm_inverse_sitediagonal_minus_input_2");
}

BOOST_AUTO_TEST_CASE( M_TM_INVERSE_SITEDIAGONAL_MINUS_3)
{
	test_m_tm_inverse_sitediagonal_minus("/m_tm_inverse_sitediagonal_minus_input_3");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(DSLASH_EO )

BOOST_AUTO_TEST_CASE( DSLASH_EO_1)
{
	test_dslash_eo("/dslash_eo_input_1");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_2)
{
	test_dslash_eo("/dslash_eo_input_2");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_3)
{
	test_dslash_eo("/dslash_eo_input_3");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_4)
{
	test_dslash_eo("/dslash_eo_input_4");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_5)
{
	test_dslash_eo("/dslash_eo_input_5");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_6)
{
	test_dslash_eo("/dslash_eo_input_6");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_7)
{
	test_dslash_eo("/dslash_eo_input_7");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_8)
{
	test_dslash_eo("/dslash_eo_input_8");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_9)
{
	test_dslash_eo("/dslash_eo_input_9");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_10)
{
	test_dslash_eo("/dslash_eo_input_10");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_11)
{
	test_dslash_eo("/dslash_eo_input_11");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_12)
{
	test_dslash_eo("/dslash_eo_input_12");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_13)
{
	test_dslash_eo("/dslash_eo_input_13");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_14)
{
	test_dslash_eo("/dslash_eo_input_14");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_15)
{
	test_dslash_eo("/dslash_eo_input_15");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_16)
{
	test_dslash_eo("/dslash_eo_input_16");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_17)
{
	test_dslash_eo("/dslash_eo_input_17");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_18)
{
	test_dslash_eo("/dslash_eo_input_18");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_19)
{
	test_dslash_eo("/dslash_eo_input_19");
}

BOOST_AUTO_TEST_CASE( DSLASH_EO_20)
{
	test_dslash_eo("/dslash_eo_input_20");
}

BOOST_AUTO_TEST_SUITE_END()

void test_m_fermion_compare_noneo_eo(std::string inputfile, int switcher)
{
	//switcher switches between similar functions
	//0: m_wilson (pure wilson)
	//1: m_tm_plus (twisted mass, upper flavour)
	//2: m_tm_minus (twisted mass, lower flavour)
	using namespace hardware::buffers;

	std::string kernelName = "Test equivalence of ";
	if(switcher == 0) {
		kernelName += "m_wilson";
	} else if(switcher == 1) {
		kernelName += "m_tm_plus";
	} else if(switcher == 2) {
		kernelName += "m_tm_minus";
	} else {
		logger.fatal() << "wrong parameter in test_m_fermion";
	}
	kernelName += " in eo- and non-eo formulation";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	cl_int err = CL_SUCCESS;
	hardware::code::Fermions * device = cpu.get_device();
	spinor * sf_in_noneo;
	spinor * sf_out_noneo;
	spinor * sf_in_eo1;
	spinor * sf_in_eo2;
	spinor * sf_out_eo;

	logger.info() << "Fill buffers...";
	size_t NUM_ELEMENTS_SF = meta::get_spinorfieldsize(params);
	size_t NUM_ELEMENTS_SF_EO = meta::get_eoprec_spinorfieldsize(params);

	sf_in_noneo = new spinor[NUM_ELEMENTS_SF];
	sf_out_noneo = new spinor[NUM_ELEMENTS_SF];
	sf_in_eo1 = new spinor[NUM_ELEMENTS_SF_EO];
	sf_in_eo2 = new spinor[NUM_ELEMENTS_SF_EO];
	sf_out_eo = new spinor[NUM_ELEMENTS_SF];

	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) {
		fill_sf_with_one(sf_in_eo1, NUM_ELEMENTS_SF_EO);
		fill_sf_with_one(sf_in_eo2, NUM_ELEMENTS_SF_EO);
		fill_sf_with_one(sf_in_noneo, NUM_ELEMENTS_SF);
	} else {
		fill_sf_with_random(sf_in_eo1, NUM_ELEMENTS_SF_EO, 123456);
		fill_sf_with_random(sf_in_eo2, NUM_ELEMENTS_SF_EO, 78910);
		fill_sf_with_random(sf_in_noneo, NUM_ELEMENTS_SF, 123456);
	}
	BOOST_REQUIRE(sf_in_eo1);
	BOOST_REQUIRE(sf_in_eo2);
	BOOST_REQUIRE(sf_in_noneo);
	fill_sf_with_one(sf_out_noneo, NUM_ELEMENTS_SF);
	fill_sf_with_one(sf_out_eo, NUM_ELEMENTS_SF);

	const Spinor in_eo1(NUM_ELEMENTS_SF_EO, device->get_device());
	const Spinor in_eo2(NUM_ELEMENTS_SF_EO, device->get_device());
	const Plain<spinor> out_eo(NUM_ELEMENTS_SF, device->get_device());
	const Plain<spinor> in_noneo(NUM_ELEMENTS_SF, device->get_device());
	const Plain<spinor> out_noneo(NUM_ELEMENTS_SF, device->get_device());
	//create 3 buffers for intermediate results
	const Spinor out_tmp_eo1(NUM_ELEMENTS_SF_EO, device->get_device());
	const Spinor out_tmp_eo2(NUM_ELEMENTS_SF_EO, device->get_device());
	const Spinor tmp_eo(NUM_ELEMENTS_SF_EO, device->get_device());
	in_eo1.load(sf_in_eo1);
	in_eo2.load(sf_in_eo2);
	out_eo.load(sf_out_eo);
	in_noneo.load(sf_in_noneo);
	out_noneo.load(sf_out_noneo);

	auto spinor_code = device->get_device()->get_spinor_code();
	auto gf_code = device->get_device()->get_gaugefield_code();

	if(params.get_solver() != meta::Inputparameters::cg) {
		//use use_pointsource to choose whether to copy the eo rnd vectors to noneo or vise versa
		if(params.get_use_pointsource())
			spinor_code->convert_from_eoprec_device(&in_eo1, &in_eo2, &in_noneo);
		else
			spinor_code->convert_to_eoprec_device(&in_eo1, &in_eo2, &in_noneo);
	}

	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	//create -1 on device
	hmc_complex minusone_tmp = { -1., 0.};
	hardware::buffers::Plain<hmc_complex> minusone(1, device->get_device());
	minusone.load(&minusone_tmp);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	logger.info() << "|phi_noneo|^2:";
	hmc_float cpu_back_noneo;
	spinor_code->set_float_to_global_squarenorm_device(&in_noneo, &sqnorm);
	sqnorm.dump(&cpu_back_noneo);
	logger.info() << cpu_back_noneo;
	logger.info() << "Run kernel";
	if(switcher == 0) {
		device->M_wilson_device(&in_noneo, &out_noneo,  gf_code->get_gaugefield(), params.get_kappa());
	} else if(switcher == 1) {
		device->M_tm_plus_device(&in_noneo, &out_noneo,  gf_code->get_gaugefield(), params.get_kappa(), meta::get_mubar(params));
	} else if(switcher == 2) {
		device->M_tm_minus_device(&in_noneo, &out_noneo,  gf_code->get_gaugefield(), params.get_kappa(), meta::get_mubar(params));
	} else {
		logger.fatal() << "wrong parameter in test_m_fermion";
	}
	logger.info() << "result:";
	hmc_float cpu_res_noneo;
	spinor_code->set_float_to_global_squarenorm_device(&out_noneo, &sqnorm);
	sqnorm.dump(&cpu_res_noneo);
	logger.info() << cpu_res_noneo;

	logger.info() << "|phi_eo1|^2:";
	hmc_float cpu_back_eo1;
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in_eo1, &sqnorm);
	sqnorm.dump(&cpu_back_eo1);
	logger.info() << cpu_back_eo1;
	logger.info() << "|phi_eo2|^2:";
	hmc_float cpu_back_eo2;
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in_eo2, &sqnorm);
	sqnorm.dump(&cpu_back_eo2);
	logger.info() << cpu_back_eo2;
	logger.info() << "Run kernel";
	if(switcher == 0) {
		//suppose in1 is the even, in2 the odd input vector
		//now calc out_tmp_eo1 = (1 in1 + D_eo in2)
		spinor_code->set_zero_spinorfield_eoprec_device(&tmp_eo);
		spinor_code->set_zero_spinorfield_eoprec_device(&out_tmp_eo1);

		device->dslash_eo_device(&in_eo2, &out_tmp_eo1, gf_code->get_gaugefield(), EO, params.get_kappa());

		spinor_code->saxpy_eoprec_device(&out_tmp_eo1, &in_eo1, &minusone, &out_tmp_eo1);

		//now calc out_tmp_eo2 = ( 1 in2 + D_oe in1)
		spinor_code->set_zero_spinorfield_eoprec_device(&tmp_eo);
		spinor_code->set_zero_spinorfield_eoprec_device(&out_tmp_eo2);

		device->dslash_eo_device(&in_eo1, &out_tmp_eo2, gf_code->get_gaugefield(), OE, params.get_kappa());

		spinor_code->saxpy_eoprec_device(&out_tmp_eo2, &in_eo2, &minusone, &out_tmp_eo2);

		//now, both output vectors have to be converted back to noneo
		spinor_code->convert_from_eoprec_device(&out_tmp_eo1, &out_tmp_eo2, &out_eo);
	} else if(switcher == 1) {
		//suppose in1 is the even, in2 the odd input vector
		//now calc out_tmp_eo1 = (R_even in1 + D_eo in2)
		spinor_code->set_zero_spinorfield_eoprec_device(&tmp_eo);
		spinor_code->set_zero_spinorfield_eoprec_device(&out_tmp_eo1);

		device->dslash_eo_device(&in_eo2, &out_tmp_eo1, gf_code->get_gaugefield(), EO, params.get_kappa());
		device->M_tm_sitediagonal_device(&in_eo1, &tmp_eo,  meta::get_mubar(params));

		spinor_code->saxpy_eoprec_device(&out_tmp_eo1, &tmp_eo, &minusone, &out_tmp_eo1);

		//now calc out_tmp_eo2 = ( R_odd in2 + D_oe in1)
		spinor_code->set_zero_spinorfield_eoprec_device(&tmp_eo);
		spinor_code->set_zero_spinorfield_eoprec_device(&out_tmp_eo2);

		device->dslash_eo_device(&in_eo1, &out_tmp_eo2, gf_code->get_gaugefield(), OE, params.get_kappa());
		device->M_tm_sitediagonal_device(&in_eo2, &tmp_eo,  meta::get_mubar(params));

		spinor_code->saxpy_eoprec_device(&out_tmp_eo2, &tmp_eo, &minusone, &out_tmp_eo2);

		//now, both output vectors have to be converted back to noneo
		spinor_code->convert_from_eoprec_device(&out_tmp_eo1, &out_tmp_eo2, &out_eo);
	} else if(switcher == 2) {
		//suppose in1 is the even, in2 the odd input vector
		//now calc out_tmp_eo1 = (R_even in1 + D_eo in2)
		spinor_code->set_zero_spinorfield_eoprec_device(&tmp_eo);
		spinor_code->set_zero_spinorfield_eoprec_device(&out_tmp_eo1);

		device->dslash_eo_device(&in_eo2, &out_tmp_eo1, gf_code->get_gaugefield(), EO, params.get_kappa());
		device->M_tm_sitediagonal_minus_device(&in_eo1, &tmp_eo,  meta::get_mubar(params));

		spinor_code->saxpy_eoprec_device(&out_tmp_eo1, &tmp_eo, &minusone, &out_tmp_eo1);

		//now calc out_tmp_eo2 = ( R_odd in2 + D_oe in1)
		spinor_code->set_zero_spinorfield_eoprec_device(&tmp_eo);
		spinor_code->set_zero_spinorfield_eoprec_device(&out_tmp_eo2);

		device->dslash_eo_device(&in_eo1, &out_tmp_eo2, gf_code->get_gaugefield(), OE, params.get_kappa());
		device->M_tm_sitediagonal_minus_device(&in_eo2, &tmp_eo,  meta::get_mubar(params));

		spinor_code->saxpy_eoprec_device(&out_tmp_eo2, &tmp_eo, &minusone, &out_tmp_eo2);

		//now, both output vectors have to be converted back to noneo
		spinor_code->convert_from_eoprec_device(&out_tmp_eo1, &out_tmp_eo2, &out_eo);
	}
	logger.info() << "result:";
	hmc_float cpu_res_eo;
	spinor_code->set_float_to_global_squarenorm_device(&out_eo, &sqnorm);
	sqnorm.dump(&cpu_res_eo);
	logger.info() << cpu_res_eo;

	logger.info() << "Finalize device";
	cpu.finalize();

	logger.info() << "Clear buffers";
	delete[] sf_in_noneo;
	delete[] sf_out_noneo;
	delete[] sf_in_eo1;
	delete[] sf_in_eo2;
	delete[] sf_out_eo;

	logger.info() << "Compare eo and non-eo results";
	BOOST_REQUIRE_CLOSE(cpu_res_eo, cpu_res_noneo, params.get_solver_prec() );
	testFloatAgainstInputparameters(cpu_res_noneo, params);
	testFloatAgainstInputparameters(cpu_res_eo, params);
	BOOST_MESSAGE("Test done");
}

void test_m_wilson_compare_noneo_eo(std::string inputfile)
{
	test_m_fermion_compare_noneo_eo(inputfile, 0);
}

void test_m_tm_plus_compare_noneo_eo(std::string inputfile)
{
	test_m_fermion_compare_noneo_eo(inputfile, 1);
}

void test_m_tm_minus_compare_noneo_eo(std::string inputfile)
{
	test_m_fermion_compare_noneo_eo(inputfile, 2);
}

BOOST_AUTO_TEST_SUITE(M_WILSON_COMPARE_NONEO_EO )

BOOST_AUTO_TEST_CASE(M_WILSON_COMPARE_NONEO_EO_1)
{
	test_m_wilson_compare_noneo_eo("/m_wilson_compare_noneo_eo_input_1");
}

BOOST_AUTO_TEST_CASE(M_WILSON_COMPARE_NONEO_EO_2)
{
	test_m_wilson_compare_noneo_eo("/m_wilson_compare_noneo_eo_input_2");
}

BOOST_AUTO_TEST_CASE(M_WILSON_COMPARE_NONEO_EO_3)
{
	test_m_wilson_compare_noneo_eo("/m_wilson_compare_noneo_eo_input_3");
}

BOOST_AUTO_TEST_CASE(M_WILSON_COMPARE_NONEO_EO_4)
{
	test_m_wilson_compare_noneo_eo("/m_wilson_compare_noneo_eo_input_4");
}

BOOST_AUTO_TEST_CASE(M_WILSON_COMPARE_NONEO_EO_5)
{
	test_m_wilson_compare_noneo_eo("/m_wilson_compare_noneo_eo_input_5");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(M_TM_COMPARE_NONEO_EO )

BOOST_AUTO_TEST_CASE(M_TM_COMPARE_NONEO_EO_1)
{
	test_m_tm_plus_compare_noneo_eo("/m_tm_compare_noneo_eo_input_1");
}

BOOST_AUTO_TEST_CASE(M_TM_COMPARE_NONEO_EO_2)
{
	test_m_tm_plus_compare_noneo_eo("/m_tm_compare_noneo_eo_input_2");
}

BOOST_AUTO_TEST_CASE(M_TM_COMPARE_NONEO_EO_3)
{
	test_m_tm_plus_compare_noneo_eo("/m_tm_compare_noneo_eo_input_3");
}

BOOST_AUTO_TEST_CASE(M_TM_COMPARE_NONEO_EO_4)
{
	test_m_tm_plus_compare_noneo_eo("/m_tm_compare_noneo_eo_input_4");
}

BOOST_AUTO_TEST_CASE(M_TM_COMPARE_NONEO_EO_5)
{
	test_m_tm_plus_compare_noneo_eo("/m_tm_compare_noneo_eo_input_5");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(M_TM_MINUS_COMPARE_NONEO_EO )

BOOST_AUTO_TEST_CASE(M_TM_MINUS_COMPARE_NONEO_EO_1)
{
	test_m_tm_plus_compare_noneo_eo("/m_tm_minus_compare_noneo_eo_input_1");
}

BOOST_AUTO_TEST_CASE(M_TM_MINUS_COMPARE_NONEO_EO_2)
{
	test_m_tm_plus_compare_noneo_eo("/m_tm_minus_compare_noneo_eo_input_2");
}

BOOST_AUTO_TEST_CASE(M_TM_MINUS_COMPARE_NONEO_EO_3)
{
	test_m_tm_plus_compare_noneo_eo("/m_tm_minus_compare_noneo_eo_input_3");
}

BOOST_AUTO_TEST_CASE(M_TM_MINUS_COMPARE_NONEO_EO_4)
{
	test_m_tm_plus_compare_noneo_eo("/m_tm_minus_compare_noneo_eo_input_4");
}

BOOST_AUTO_TEST_CASE(M_TM_MINUS_COMPARE_NONEO_EO_5)
{
	test_m_tm_plus_compare_noneo_eo("/m_tm_minus_compare_noneo_eo_input_5");
}

BOOST_AUTO_TEST_SUITE_END()
