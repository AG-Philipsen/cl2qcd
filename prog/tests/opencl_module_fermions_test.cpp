#include "../opencl_module_fermions.h"
#include "../gaugefield_hybrid.h"
#include "../meta/util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE_FERMIONS
#include <boost/test/unit_test.hpp>

//some functionality
#include "test_util.h"

extern std::string const version;
std::string const version = "0.1";

class TestGaugefield : public Gaugefield_hybrid {

public:
	TestGaugefield(const hardware::System * system) : Gaugefield_hybrid(system) {
		auto inputfile = system->get_inputparameters();
		init(1, inputfile.get_use_gpu() ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU);
		meta::print_info_hmc("test program", inputfile);
	};

	virtual void init_tasks();
	virtual void finalize_opencl();

	Opencl_Module_Fermions * get_device();
};

void TestGaugefield::init_tasks()
{
	opencl_modules = new Opencl_Module* [get_num_tasks()];
	opencl_modules[0] = new Opencl_Module_Fermions(get_parameters(), get_device_for_task(0));
	opencl_modules[0]->init();
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

void fill_sf_with_random(spinor * sf_in, int size)
{
	prng_init(123456);
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

Opencl_Module_Fermions* TestGaugefield::get_device()
{
	return static_cast<Opencl_Module_Fermions*>(opencl_modules[0]);
}

void test_m_tm_plus(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName = "m_tm_plus";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	cl_int err = CL_SUCCESS;
	Opencl_Module_Fermions * device = cpu.get_device();
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

	logger.info() << "|phi|^2:";
	hmc_float cpu_back;
	device->set_float_to_global_squarenorm_device(&in, &sqnorm);
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;
	logger.info() << "Run kernel";
	device->M_tm_plus_device(&in, &out,  device->get_gaugefield(), params.get_kappa(), meta::get_mubar(params));
	logger.info() << "result:";
	hmc_float cpu_res;
	device->set_float_to_global_squarenorm_device(&out, &sqnorm);
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

void test_dslash_eo(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName = "dslash_eo";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	Opencl_Module_Fermions * device = cpu.get_device();
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
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	size_t NUM_ELEMENTS_SF_EO = meta::get_eoprec_spinorfieldsize(params);
	const Spinor in_eo_even(NUM_ELEMENTS_SF_EO, device->get_device());
	const Spinor in_eo_odd(NUM_ELEMENTS_SF_EO, device->get_device());
	const Spinor out_eo(NUM_ELEMENTS_SF_EO, device->get_device());
	device->convert_to_eoprec_device(&in_eo_even, &in_eo_odd, &in);

	logger.info() << "|phi|^2:";
	hmc_float cpu_back;
	device->set_float_to_global_squarenorm_eoprec_device(&in_eo_even, &sqnorm);
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;

	//switch according to "use_pointsource"
	hmc_float cpu_res;
	if(params.get_use_pointsource()) {
		device->dslash_eo_device( &in_eo_even, &out_eo, device->get_gaugefield(), EVEN, params.get_kappa() );
	} else {
		device->dslash_eo_device( &in_eo_even, &out_eo, device->get_gaugefield(), ODD, params.get_kappa() );
	}
	device->set_float_to_global_squarenorm_eoprec_device(&out_eo, &sqnorm);
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


BOOST_AUTO_TEST_SUITE( M_WILSON )

BOOST_AUTO_TEST_CASE( M_WILSON_1)
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( M_TM_MINUS  )

BOOST_AUTO_TEST_CASE( M_TM_MINUS_1 )
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!");
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
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( GAMMA5 )

BOOST_AUTO_TEST_CASE( GAMMA5_1)
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( GAMMA5_EO)

BOOST_AUTO_TEST_CASE( GAMMA5_EO_1)
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(M_TM_SITEDIAGONAL )

BOOST_AUTO_TEST_CASE( M_TM_SITEDIAGONAL1_)
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( M_TM_INVERSE_SITEDIAGONAL )

BOOST_AUTO_TEST_CASE( M_TM_INVERSE_SITEDIAGONAL_1)
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(M_TM_SITEDIAGONAL_MINUS )

BOOST_AUTO_TEST_CASE( M_TM_SITEDIAGONAL_MINUS_1)
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( M_TM_INVERSE_SITEDIAGONAL_MINUS)

BOOST_AUTO_TEST_CASE( M_TM_INVERSE_SITEDIAGONAL_MINUS_1 )
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!");
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

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(DSLASH_AND_GAMMA5_EO )

BOOST_AUTO_TEST_CASE( DSLASH_AND_GAMMA5_EO_1)
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_EO )

BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_EO_1)
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_MINUS_EO )

BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_MINUS_EO_1)
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(M_TM_SITEDIAGONAL_AND_GAMMA5_EO )

BOOST_AUTO_TEST_CASE(M_TM_SITEDIAGONAL_AND_GAMMA5_EO_1)
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(M_TM_SITEDIAGONAL_MINUS_AND_GAMMA5_EO )

BOOST_AUTO_TEST_CASE(M_TM_SITEDIAGONAL_MINUS_AND_GAMMA5_EO_1)
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

