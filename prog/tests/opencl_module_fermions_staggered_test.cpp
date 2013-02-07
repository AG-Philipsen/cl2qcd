#include "../meta/util.hpp"
#include "../host_random.h"
#include "../physics/lattices/gaugefield.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE_FERMIONS_STAGGERED
#include <boost/test/unit_test.hpp>

//some functionality
#include "test_util.h"

class TestGaugefield {

public:
	TestGaugefield(const hardware::System * system) : system(system), prng(*system), gf(*system, prng) {
		BOOST_REQUIRE_EQUAL(system->get_devices().size(), 1);
		auto inputfile = system->get_inputparameters();
		meta::print_info_hmc("test program", inputfile);
	};

	hardware::code::Fermions_staggered * get_device();
	const hardware::buffers::SU3 * get_gaugefield();

private:
	const hardware::System * const system;
	physics::PRNG prng;
	const physics::lattices::Gaugefield gf;
};

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

hardware::code::Fermions_staggered* TestGaugefield::get_device()
{
	return system->get_devices()[0]->get_fermion_staggered_code();
}

const hardware::buffers::SU3 * TestGaugefield::get_gaugefield()
{
	return gf.get_buffers().at(0);
}

void test_build(std::string inputfile)
{
	logger.info() << "build opencl_module_fermions_staggered";
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	BOOST_MESSAGE("Test done");
}

void test_m_staggered(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName;
	kernelName = "m_staggered";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	cl_int err = CL_SUCCESS;
	hardware::code::Fermions_staggered * device = cpu.get_device();
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

	logger.info() << "|phi|^2:";
	hmc_float cpu_back;
	spinor_code->set_float_to_global_squarenorm_device(&in, &sqnorm);
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;
	logger.info() << "Run kernel";
	device->M_staggered_device(&in, &out,  cpu.get_gaugefield(), params.get_kappa());
	logger.info() << "result:";
	hmc_float cpu_res;
	spinor_code->set_float_to_global_squarenorm_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	logger.info() << "Clear buffers";
	delete[] sf_in;
	delete[] sf_out;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

BOOST_AUTO_TEST_SUITE(BUILD)

BOOST_AUTO_TEST_CASE( BUILD_1 )
{
	test_build("/opencl_module_fermions_staggered_build_input_1");
}

BOOST_AUTO_TEST_CASE( BUILD_2 )
{
	test_build("/opencl_module_fermions_staggered_build_input_2");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( M_STAGGERED )

BOOST_AUTO_TEST_CASE( M_STAGGERED_1)
{
	test_m_staggered("/m_staggered_input_1");
}
/**
 * to be added...
 */
/*
BOOST_AUTO_TEST_CASE( M_STAGGERED_2)
{
	test_m_staggered("/m_staggered_input_2");
}

BOOST_AUTO_TEST_CASE( M_STAGGERED_3)
{
	test_m_staggered("/m_staggered_input_3");
}

BOOST_AUTO_TEST_CASE( M_STAGGERED_4)
{
	test_m_staggered("/m_staggered_input_4");
}

BOOST_AUTO_TEST_CASE( M_STAGGERED_5)
{
	test_m_staggered("/m_staggered_input_5");
}

BOOST_AUTO_TEST_CASE( M_STAGGERED_6)
{
	test_m_staggered("/m_staggered_input_6");
}
*/
BOOST_AUTO_TEST_SUITE_END()

