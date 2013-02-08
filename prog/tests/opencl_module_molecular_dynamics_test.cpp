#include "../meta/util.hpp"
#include "../physics/lattices/gaugefield.hpp"
#include "../hardware/device.hpp"
#include "../hardware/code/molecular_dynamics.hpp"
#include "../hardware/code/spinors.hpp"
#include "../hardware/code/gaugemomentum.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE_MOLECULAR_DYNAMICS
#include <boost/test/unit_test.hpp>

#include "test_util.h"
#include "../host_random.h"

class TestGaugefield {

public:
	TestGaugefield(const hardware::System * system) : system(system), prng(*system), gf(*system, prng) {
		BOOST_REQUIRE_EQUAL(system->get_devices().size(), 1);
		auto inputfile = system->get_inputparameters();
		meta::print_info_hmc("test program", inputfile);
	};

	const hardware::code::Molecular_Dynamics * get_device();
	const hardware::buffers::SU3 * get_gaugefield();

	void print_gaugeobservables() {
		physics::lattices::print_gaugeobservables(gf, 0);
	}

private:
	const hardware::System * const system;
	physics::PRNG prng;
	const physics::lattices::Gaugefield gf;
};

const hardware::code::Molecular_Dynamics* TestGaugefield::get_device()
{
	return system->get_devices()[0]->get_molecular_dynamics_code();
}

const hardware::buffers::SU3 * TestGaugefield::get_gaugefield()
{
	return gf.get_buffers().at(0);
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

void fill_with_one(hmc_float * sf_in, int size)
{
	for(int i = 0; i < size; ++i) {
		sf_in[i] = 1.;
	}
	return;
}

void fill_with_random(hmc_float * sf_in, int size, int seed)
{
	prng_init(seed);
	for(int i = 0; i < size; ++i) {
		sf_in[i] = prng_double();
	}
	return;
}

ae make_ae(hmc_float e1, hmc_float e2, hmc_float e3, hmc_float e4,
           hmc_float e5, hmc_float e6, hmc_float e7, hmc_float e8)
{
	ae tmp = {e1, e2, e3, e4, e5, e6, e7, e8};
	return tmp;
}

void fill_with_zero(ae * ae, int size)
{
	for(int i = 0; i < size; ++i) {
		ae[i] = make_ae(0., 0., 0., 0., 0., 0., 0., 0.);
	}
	return;
}

void fill_sf_with_random(spinor * sf_in, int size, int seed)
{
	//  Random rnd_loc(seed);
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

void fill_sf_with_random_eo(spinor * sf_in1, spinor * sf_in2, int size, int seed)
{
	prng_init(seed);
	for(int i = 0; i < size; ++i) {
		sf_in1[i].e0.e0.re = prng_double();
		sf_in1[i].e0.e1.re = prng_double();
		sf_in1[i].e0.e2.re = prng_double();
		sf_in1[i].e1.e0.re = prng_double();
		sf_in1[i].e1.e1.re = prng_double();
		sf_in1[i].e1.e2.re = prng_double();
		sf_in1[i].e2.e0.re = prng_double();
		sf_in1[i].e2.e1.re = prng_double();
		sf_in1[i].e2.e2.re = prng_double();
		sf_in1[i].e3.e0.re = prng_double();
		sf_in1[i].e3.e1.re = prng_double();
		sf_in1[i].e3.e2.re = prng_double();

		sf_in1[i].e0.e0.im = prng_double();
		sf_in1[i].e0.e1.im = prng_double();
		sf_in1[i].e0.e2.im = prng_double();
		sf_in1[i].e1.e0.im = prng_double();
		sf_in1[i].e1.e1.im = prng_double();
		sf_in1[i].e1.e2.im = prng_double();
		sf_in1[i].e2.e0.im = prng_double();
		sf_in1[i].e2.e1.im = prng_double();
		sf_in1[i].e2.e2.im = prng_double();
		sf_in1[i].e3.e0.im = prng_double();
		sf_in1[i].e3.e1.im = prng_double();
		sf_in1[i].e3.e2.im = prng_double();

		sf_in2[i].e0.e0.re = prng_double();
		sf_in2[i].e0.e1.re = prng_double();
		sf_in2[i].e0.e2.re = prng_double();
		sf_in2[i].e1.e0.re = prng_double();
		sf_in2[i].e1.e1.re = prng_double();
		sf_in2[i].e1.e2.re = prng_double();
		sf_in2[i].e2.e0.re = prng_double();
		sf_in2[i].e2.e1.re = prng_double();
		sf_in2[i].e2.e2.re = prng_double();
		sf_in2[i].e3.e0.re = prng_double();
		sf_in2[i].e3.e1.re = prng_double();
		sf_in2[i].e3.e2.re = prng_double();

		sf_in2[i].e0.e0.im = prng_double();
		sf_in2[i].e0.e1.im = prng_double();
		sf_in2[i].e0.e2.im = prng_double();
		sf_in2[i].e1.e0.im = prng_double();
		sf_in2[i].e1.e1.im = prng_double();
		sf_in2[i].e1.e2.im = prng_double();
		sf_in2[i].e2.e0.im = prng_double();
		sf_in2[i].e2.e1.im = prng_double();
		sf_in2[i].e2.e2.im = prng_double();
		sf_in2[i].e3.e0.im = prng_double();
		sf_in2[i].e3.e1.im = prng_double();
		sf_in2[i].e3.e2.im = prng_double();
	}
	return;
}

void test_build(std::string inputfile)
{
	logger.info() << "build opencl_module_molecular_dynamics";
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	BOOST_MESSAGE("Test done");
}

void test_stout_smear_fermion_force(std::string inputfile)
{

}

void test_gf_update(std::string inputfile)
{
	std::string kernelName = "md_update_gaugefield";
	printKernelInfo(kernelName);

	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	auto * device = cpu.get_device();
	hmc_float * gm_in;

	logger.info() << "create buffers";
	size_t NUM_ELEMENTS_AE = meta::get_vol4d(params) * NDIM * meta::get_su3algebrasize();
	gm_in = new hmc_float[NUM_ELEMENTS_AE];

	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) {
		fill_with_one(gm_in, NUM_ELEMENTS_AE);
	} else {
		fill_with_random(gm_in, NUM_ELEMENTS_AE, 123456);
	}
	BOOST_REQUIRE(gm_in);
	auto gf_code = device->get_device()->get_gaugefield_code();
	auto gm_code = device->get_device()->get_gaugemomentum_code();

	hardware::buffers::Gaugemomentum in(meta::get_vol4d(params) * NDIM, device->get_device());
	gm_code->importGaugemomentumBuffer(&in, reinterpret_cast<ae*>(gm_in));
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	logger.info() << "|in|^2:";

	gm_code->set_float_to_gaugemomentum_squarenorm_device(&in, &sqnorm);
	hmc_float cpu_back;
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;

	logger.info() << "Run kernel";
	hmc_float eps = params.get_tau();
	device->md_update_gaugefield_device(&in, cpu.get_gaugefield(), eps);
	logger.info() << "gaugeobservables: ";
	cpu.print_gaugeobservables();

	hmc_float plaq_cpu, tplaq_cpu, splaq_cpu;
	hmc_complex pol_cpu;
	gf_code->gaugeobservables(cpu.get_gaugefield(), &plaq_cpu, &tplaq_cpu, &splaq_cpu, &pol_cpu);

	logger.info() << "Free buffers";
	delete[] gm_in;

	testFloatAgainstInputparameters(plaq_cpu, params);
	BOOST_MESSAGE("Test done");
}

void test_f_update(std::string inputfile)
{
	std::string kernelName = "md_update_gaugemomenta";
	printKernelInfo(kernelName);

	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	auto * device = cpu.get_device();
	hmc_float * gm_in;
	hmc_float * gm_out;
	auto gm_code = device->get_device()->get_gaugemomentum_code();

	logger.info() << "create buffers";
	size_t NUM_ELEMENTS_AE = meta::get_vol4d(params) * NDIM * meta::get_su3algebrasize();
	gm_in = new hmc_float[NUM_ELEMENTS_AE];
	gm_out = new hmc_float[NUM_ELEMENTS_AE];

	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) {
		fill_with_one(gm_in, NUM_ELEMENTS_AE);
		fill_with_one(gm_out, NUM_ELEMENTS_AE);
	} else {
		fill_with_random(gm_in, NUM_ELEMENTS_AE, 123456);
		fill_with_random(gm_out, NUM_ELEMENTS_AE, 789101);
	}
	BOOST_REQUIRE(gm_in);
	BOOST_REQUIRE(gm_out);

	hardware::buffers::Gaugemomentum in(meta::get_vol4d(params) * NDIM, device->get_device());
	gm_code->importGaugemomentumBuffer(&in, reinterpret_cast<ae*>(gm_in));
	hardware::buffers::Gaugemomentum out(meta::get_vol4d(params) * NDIM, device->get_device());
	gm_code->importGaugemomentumBuffer(&out, reinterpret_cast<ae*>(gm_out));
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	logger.info() << "|in|^2:";
	gm_code->set_float_to_gaugemomentum_squarenorm_device(&in, &sqnorm);
	hmc_float cpu_back;
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;
	logger.info() << "|out|^2:";
	gm_code->set_float_to_gaugemomentum_squarenorm_device(&out, &sqnorm);
	hmc_float cpu_back2;
	sqnorm.dump(&cpu_back2);
	logger.info() << cpu_back2;

	logger.info() << "Run kernel";
	hmc_float eps = params.get_tau();
	device->md_update_gaugemomentum_device(&in, &out, eps);

	gm_code->set_float_to_gaugemomentum_squarenorm_device(&out, &sqnorm);
	hmc_float cpu_res;
	sqnorm.dump(&cpu_res);
	logger.info() << "result:";
	logger.info() << cpu_res;

	logger.info() << "Free buffers";
	delete[] gm_in;
	delete[] gm_out;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_f_gauge(std::string inputfile)
{
	std::string kernelName = "f_gauge";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	auto * device = cpu.get_device();
	ae * gm_out;
	auto gm_code = device->get_device()->get_gaugemomentum_code();

	logger.info() << "create buffers";
	size_t NUM_ELEMENTS_AE = meta::get_vol4d(params) * NDIM * meta::get_su3algebrasize();
	gm_out = new ae[NUM_ELEMENTS_AE];
	fill_with_zero(gm_out, NUM_ELEMENTS_AE);
	BOOST_REQUIRE(gm_out);

	hardware::buffers::Gaugemomentum out(meta::get_vol4d(params) * NDIM, device->get_device());
	gm_code->importGaugemomentumBuffer(&out, reinterpret_cast<ae*>(gm_out));
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	logger.info() << "|out|^2:";
	gm_code->set_float_to_gaugemomentum_squarenorm_device(&out, &sqnorm);
	hmc_float cpu_back;
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;

	device->gauge_force_device( cpu.get_gaugefield(), &out);

	logger.info() << "result:";
	hmc_float cpu_res;
	gm_code->set_float_to_gaugemomentum_squarenorm_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");

	logger.info() << "Finalize device";
	delete[] gm_out;
}

void test_f_gauge_tlsym(std::string inputfile)
{
	std::string kernelName = "f_gauge_tlsym";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	auto * device = cpu.get_device();
	ae * gm_out;
	auto gm_code = device->get_device()->get_gaugemomentum_code();

	logger.info() << "create buffers";
	size_t NUM_ELEMENTS_AE = meta::get_vol4d(params) * NDIM * meta::get_su3algebrasize();
	gm_out = new ae[NUM_ELEMENTS_AE];
	fill_with_zero(gm_out, NUM_ELEMENTS_AE);
	BOOST_REQUIRE(gm_out);

	hardware::buffers::Gaugemomentum out(meta::get_vol4d(params) * NDIM, device->get_device());
	gm_code->importGaugemomentumBuffer(&out, reinterpret_cast<ae*>(gm_out));
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	logger.info() << "|out|^2:";
	gm_code->set_float_to_gaugemomentum_squarenorm_device(&out, &sqnorm);
	hmc_float cpu_back;
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;

	device->gauge_force_tlsym_device( cpu.get_gaugefield(), &out);

	logger.info() << "result:";
	hmc_float cpu_res;
	gm_code->set_float_to_gaugemomentum_squarenorm_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");

	delete[] gm_out;
}

void test_f_fermion(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName = "f_fermion";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	auto * device = cpu.get_device();
	spinor * sf_in1;
	spinor * sf_in2;
	ae * ae_out;

	logger.info() << "Fill buffers...";
	int NUM_ELEMENTS_SF = meta::get_spinorfieldsize(params);
	int NUM_ELEMENTS_AE = meta::get_vol4d(params) * NDIM * meta::get_su3algebrasize();

	sf_in1 = new spinor[NUM_ELEMENTS_SF];
	sf_in2 = new spinor[NUM_ELEMENTS_SF];
	ae_out = new ae[NUM_ELEMENTS_AE];

	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) {
		fill_sf_with_one(sf_in1, NUM_ELEMENTS_SF);
		fill_sf_with_one(sf_in2, NUM_ELEMENTS_SF);
	} else {
		fill_sf_with_random(sf_in1, NUM_ELEMENTS_SF, 123456);
		fill_sf_with_random(sf_in2, NUM_ELEMENTS_SF, 789101);
	}
	fill_with_zero(ae_out, NUM_ELEMENTS_AE);
	BOOST_REQUIRE(sf_in1);
	BOOST_REQUIRE(sf_in2);
	BOOST_REQUIRE(ae_out);

	auto spinor_code = device->get_device()->get_spinor_code();
	auto gm_code = device->get_device()->get_gaugemomentum_code();

	const Plain<spinor> in1(NUM_ELEMENTS_SF, device->get_device());
	const Plain<spinor> in2(NUM_ELEMENTS_SF, device->get_device());
	in1.load(sf_in1);
	in2.load(sf_in2);
	Gaugemomentum out(meta::get_vol4d(params) * NDIM, device->get_device());
	gm_code->importGaugemomentumBuffer(&out, reinterpret_cast<ae*>(ae_out));
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	logger.info() << "|phi_1|^2:";
	hmc_float cpu_back;
	spinor_code->set_float_to_global_squarenorm_device(&in1, &sqnorm);
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;
	logger.info() << "|phi_2|^2:";
	hmc_float cpu_back2;
	spinor_code->set_float_to_global_squarenorm_device(&in2, &sqnorm);
	sqnorm.dump(&cpu_back2);
	logger.info() << cpu_back2;
	logger.info() << "Run kernel";
	device->fermion_force_device( &in1, &in2, cpu.get_gaugefield(), &out, params.get_kappa());
	logger.info() << "result:";
	hmc_float cpu_res;
	gm_code->set_float_to_gaugemomentum_squarenorm_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	logger.info() << "Clear buffers";
	delete[] sf_in1;
	delete[] sf_in2;
	delete[] ae_out;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

void test_f_fermion_eo(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName = "f_fermion_eo";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	auto * device = cpu.get_device();
	spinor * sf_in1;
	spinor * sf_in2;
	ae * ae_out;

	logger.info() << "fill buffers";
	size_t NUM_ELEMENTS_SF =  meta::get_eoprec_spinorfieldsize(params);
	size_t NUM_ELEMENTS_AE = meta::get_vol4d(params) * NDIM * meta::get_su3algebrasize();

	sf_in1 = new spinor[NUM_ELEMENTS_SF];
	sf_in2 = new spinor[NUM_ELEMENTS_SF];
	ae_out = new ae[NUM_ELEMENTS_AE];

	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) {
		fill_sf_with_one(sf_in1, NUM_ELEMENTS_SF);
		fill_sf_with_one(sf_in2, NUM_ELEMENTS_SF);
	} else {
		fill_sf_with_random_eo(sf_in1, sf_in2, NUM_ELEMENTS_SF, 123456);
	}
	fill_with_zero(ae_out, NUM_ELEMENTS_AE);
	BOOST_REQUIRE(sf_in1);
	BOOST_REQUIRE(sf_in2);
	BOOST_REQUIRE(ae_out);

	auto spinor_code = device->get_device()->get_spinor_code();
	auto gm_code = device->get_device()->get_gaugemomentum_code();

	const Spinor in1(NUM_ELEMENTS_SF, device->get_device());
	const Spinor in2(NUM_ELEMENTS_SF, device->get_device());
	Gaugemomentum out(meta::get_vol4d(params) * NDIM, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	spinor_code->copy_to_eoprec_spinorfield_buffer(&in1, sf_in1);
	spinor_code->copy_to_eoprec_spinorfield_buffer(&in2, sf_in2);
	gm_code->importGaugemomentumBuffer(&out, ae_out);


	hmc_float cpu_res, cpu_back, cpu_back2;
	logger.info() << "|in_1|^2:";
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in1, &sqnorm);
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;
	logger.info() << "|in_2|^2:";
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in2, &sqnorm);
	sqnorm.dump(&cpu_back2);
	logger.info() << cpu_back2;
	logger.info() << "Run kernel";

	//switch according to "read_multiple_configs"
	if(params.get_read_multiple_configs()) {
		int tmp = EVEN;
		device->fermion_force_eo_device(&in1, &in2, cpu.get_gaugefield(), &out, tmp, params.get_kappa() );
	} else {
		int tmp = ODD;
		device->fermion_force_eo_device(&in1, &in2, cpu.get_gaugefield(), &out, tmp, params.get_kappa() );
	}
	logger.info() << "|force|^2:";
	gm_code->set_float_to_gaugemomentum_squarenorm_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	logger.info() << "Clear buffers";
	delete[] sf_in1;
	delete[] sf_in2;
	delete[] ae_out;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

BOOST_AUTO_TEST_SUITE(BUILD)

BOOST_AUTO_TEST_CASE( BUILD_1 )
{
	test_build("/opencl_module_molecular_dynamics_build_input_1");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( STOUT_SMEAR_FERMION_FORCE )

BOOST_AUTO_TEST_CASE(STOUT_SMEAR_FERMION_FORCE_1 )
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!!");
	test_stout_smear_fermion_force("/stout_smear_fermion_force_input_1");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( GF_UPDATE )

BOOST_AUTO_TEST_CASE(GF_UPDATE_1 )
{
	test_gf_update("/gf_update_input_1");
}

BOOST_AUTO_TEST_CASE( GF_UPDATE_2 )
{
	test_gf_update("/gf_update_input_2");
}

BOOST_AUTO_TEST_CASE( GF_UPDATE_3 )
{
	test_gf_update("/gf_update_input_3");
}

BOOST_AUTO_TEST_CASE( GF_UPDATE_4 )
{
	test_gf_update("/gf_update_input_4");
}

BOOST_AUTO_TEST_CASE(GF_UPDATE_5 )
{
	BOOST_MESSAGE("THIS TEST HAS TO BE INVESTIGATED!! IT COULD HINT TO AN ERROR IN THE FUNCTION!");
	test_gf_update("/gf_update_input_5");
}

BOOST_AUTO_TEST_CASE(GF_UPDATE_6 )
{
	test_gf_update("/gf_update_input_6");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_UPDATE )

BOOST_AUTO_TEST_CASE( F_UPDATE_1 )
{
	test_f_update("/f_update_input_1");
}

BOOST_AUTO_TEST_CASE( F_UPDATE_2 )
{
	test_f_update("/f_update_input_2");
}

BOOST_AUTO_TEST_CASE( F_UPDATE_3 )
{
	test_f_update("/f_update_input_3");
}

BOOST_AUTO_TEST_CASE( F_UPDATE_4 )
{
	test_f_update("/f_update_input_4");
}

BOOST_AUTO_TEST_CASE( F_UPDATE_5 )
{
	test_f_update("/f_update_input_5");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_GAUGE )

BOOST_AUTO_TEST_CASE( F_GAUGE_1 )
{
	test_f_gauge("/f_gauge_input_1");
}

BOOST_AUTO_TEST_CASE( F_GAUGE_2 )
{
	test_f_gauge("/f_gauge_input_2");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_GAUGE_TLSYM )

BOOST_AUTO_TEST_CASE( F_GAUGE_TLSYM_1 )
{
	test_f_gauge_tlsym("/f_gauge_tlsym_input_1");
}

BOOST_AUTO_TEST_CASE( F_GAUGE_TLSYM_2 )
{
	test_f_gauge_tlsym("/f_gauge_tlsym_input_2");
}

BOOST_AUTO_TEST_CASE( F_GAUGE_TLSYM_3 )
{
	test_f_gauge_tlsym("/f_gauge_tlsym_input_3");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_FERMION )

BOOST_AUTO_TEST_CASE( F_FERMION_1 )
{
	test_f_fermion("/f_fermion_input_1");
}

BOOST_AUTO_TEST_CASE( F_FERMION_2 )
{
	test_f_fermion("/f_fermion_input_2");
}

BOOST_AUTO_TEST_CASE( F_FERMION_3 )
{
	test_f_fermion("/f_fermion_input_3");
}

BOOST_AUTO_TEST_CASE( F_FERMION_4 )
{
	test_f_fermion("/f_fermion_input_4");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_FERMION_EO )

BOOST_AUTO_TEST_CASE( F_FERMION_EO_1 )
{
	test_f_fermion_eo("/f_fermion_eo_input_1");
}

BOOST_AUTO_TEST_CASE( F_FERMION_EO_2 )
{
	test_f_fermion_eo("/f_fermion_eo_input_2");
}

BOOST_AUTO_TEST_CASE( F_FERMION_EO_3 )
{
	test_f_fermion_eo("/f_fermion_eo_input_3");
}

BOOST_AUTO_TEST_CASE( F_FERMION_EO_4 )
{
	test_f_fermion_eo("/f_fermion_eo_input_4");
}

BOOST_AUTO_TEST_CASE( F_FERMION_EO_5 )
{
	test_f_fermion_eo("/f_fermion_eo_input_5");
}

BOOST_AUTO_TEST_CASE( F_FERMION_EO_6 )
{
	test_f_fermion_eo("/f_fermion_eo_input_6");
}

BOOST_AUTO_TEST_CASE( F_FERMION_EO_7 )
{
	test_f_fermion_eo("/f_fermion_eo_input_7");
}

BOOST_AUTO_TEST_CASE( F_FERMION_EO_8 )
{
	test_f_fermion_eo("/f_fermion_eo_input_8");
}

BOOST_AUTO_TEST_CASE( F_FERMION_EO_9 )
{
	test_f_fermion_eo("/f_fermion_eo_input_9");
}

BOOST_AUTO_TEST_CASE( F_FERMION_EO_10 )
{
	test_f_fermion_eo("/f_fermion_eo_input_10");
}

BOOST_AUTO_TEST_CASE( F_FERMION_EO_11)
{
	test_f_fermion_eo("/f_fermion_eo_input_11");
}

BOOST_AUTO_TEST_CASE( F_FERMION_EO_12 )
{
	test_f_fermion_eo("/f_fermion_eo_input_12");
}

BOOST_AUTO_TEST_CASE( F_FERMION_EO_13)
{
	test_f_fermion_eo("/f_fermion_eo_input_13");
}

BOOST_AUTO_TEST_CASE( F_FERMION_EO_14 )
{
	test_f_fermion_eo("/f_fermion_eo_input_14");
}

BOOST_AUTO_TEST_CASE( F_FERMION_EO_15 )
{
	test_f_fermion_eo("/f_fermion_eo_input_15");
}

BOOST_AUTO_TEST_CASE( F_FERMION_EO_16 )
{
	test_f_fermion_eo("/f_fermion_eo_input_16");
}

BOOST_AUTO_TEST_CASE( F_FERMION_EO_17 )
{
	test_f_fermion_eo("/f_fermion_eo_input_17");
}

BOOST_AUTO_TEST_CASE( F_FERMION_EO_18 )
{
	test_f_fermion_eo("/f_fermion_eo_input_18");
}

BOOST_AUTO_TEST_CASE( F_FERMION_EO_19 )
{
	test_f_fermion_eo("/f_fermion_eo_input_19");
}

BOOST_AUTO_TEST_CASE( F_FERMION_EO_20 )
{
	test_f_fermion_eo("/f_fermion_eo_input_20");
}

BOOST_AUTO_TEST_SUITE_END()

void test_f_fermion_compare_noneo_eo(std::string inputfile)
{
	using namespace hardware::buffers;

	std::string kernelName = "Compare f_fermion in eo and noneo version";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	auto * device = cpu.get_device();

	spinor * sf_in1_noneo;
	spinor * sf_in2_noneo;
	ae * sf_out_noneo;
	ae * sf_out_eo;
	spinor * sf_in1_eo;
	spinor * sf_in2_eo;
	spinor * sf_in3_eo;
	spinor * sf_in4_eo;

	size_t NUM_ELEMENTS_SF_EO =  meta::get_eoprec_spinorfieldsize(params);
	size_t NUM_ELEMENTS_SF =  meta::get_spinorfieldsize(params);
	size_t NUM_ELEMENTS_AE = meta::get_vol4d(params) * NDIM;

	sf_in1_eo = new spinor[NUM_ELEMENTS_SF_EO];
	sf_in2_eo = new spinor[NUM_ELEMENTS_SF_EO];
	sf_in3_eo = new spinor[NUM_ELEMENTS_SF_EO];
	sf_in4_eo = new spinor[NUM_ELEMENTS_SF_EO];
	sf_out_eo = new ae[NUM_ELEMENTS_AE];
	sf_in1_noneo = new spinor[NUM_ELEMENTS_SF];
	sf_in2_noneo = new spinor[NUM_ELEMENTS_SF];
	sf_out_noneo = new ae[NUM_ELEMENTS_AE];

	fill_with_zero(sf_out_eo, NUM_ELEMENTS_AE);
	fill_with_zero(sf_out_noneo, NUM_ELEMENTS_AE);
	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) {
		fill_sf_with_one(sf_in1_eo, NUM_ELEMENTS_SF_EO);
		fill_sf_with_one(sf_in2_eo, NUM_ELEMENTS_SF_EO);
		fill_sf_with_one(sf_in3_eo, NUM_ELEMENTS_SF_EO);
		fill_sf_with_one(sf_in4_eo, NUM_ELEMENTS_SF_EO);
		fill_sf_with_one(sf_in1_noneo, NUM_ELEMENTS_SF);
		fill_sf_with_one(sf_in2_noneo, NUM_ELEMENTS_SF);
	} else {
		fill_sf_with_random_eo(sf_in1_eo, sf_in2_eo, NUM_ELEMENTS_SF_EO, 123456);
		fill_sf_with_random_eo(sf_in3_eo, sf_in4_eo, NUM_ELEMENTS_SF_EO, 789101);
		fill_sf_with_random(sf_in1_noneo, NUM_ELEMENTS_SF, 123456);
		fill_sf_with_random(sf_in2_noneo, NUM_ELEMENTS_SF, 789101);
	}

	BOOST_REQUIRE(sf_in1_eo);
	BOOST_REQUIRE(sf_in2_eo);
	BOOST_REQUIRE(sf_in3_eo);
	BOOST_REQUIRE(sf_in4_eo);
	BOOST_REQUIRE(sf_in1_noneo);
	BOOST_REQUIRE(sf_in2_noneo);
	BOOST_REQUIRE(sf_out_noneo);
	BOOST_REQUIRE(sf_out_eo);

	auto spinor_code = device->get_device()->get_spinor_code();
	auto gm_code = device->get_device()->get_gaugemomentum_code();

	const Spinor in1_eo(NUM_ELEMENTS_SF_EO, device->get_device());
	const Spinor in2_eo(NUM_ELEMENTS_SF_EO, device->get_device());
	const Spinor in3_eo(NUM_ELEMENTS_SF_EO, device->get_device());
	const Spinor in4_eo(NUM_ELEMENTS_SF_EO, device->get_device());
	const Plain<spinor> in1_noneo(NUM_ELEMENTS_SF, device->get_device());
	const Plain<spinor> in2_noneo(NUM_ELEMENTS_SF, device->get_device());
	const Gaugemomentum out_noneo(NUM_ELEMENTS_AE, device->get_device());
	const Gaugemomentum out_eo(NUM_ELEMENTS_AE, device->get_device());
	const Plain<hmc_float> sqnorm(1, device->get_device());

	gm_code->importGaugemomentumBuffer(&out_eo, sf_out_eo);
	gm_code->importGaugemomentumBuffer(&out_noneo, sf_out_noneo);

	//in case of rnd input, it is nontrivial to supply the same rnd vectors as eo and noneo input.
	//therefore, simply convert the eo input back to noneo
	if(params.get_solver() == meta::Inputparameters::cg) {
		spinor_code->copy_to_eoprec_spinorfield_buffer(&in1_eo, sf_in1_eo);
		spinor_code->copy_to_eoprec_spinorfield_buffer(&in2_eo, sf_in2_eo);
		spinor_code->copy_to_eoprec_spinorfield_buffer(&in3_eo, sf_in3_eo);
		spinor_code->copy_to_eoprec_spinorfield_buffer(&in4_eo, sf_in4_eo);
		in1_noneo.load(sf_in1_noneo);
		in2_noneo.load(sf_in2_noneo);
	} else {
		//one can either convert to or from eoprec, use read_multiple_configs for that
		//NOTE: there is machinery to compare vectors in the old executable
		if(params.get_read_multiple_configs()) {
			spinor_code->copy_to_eoprec_spinorfield_buffer(&in1_eo, sf_in1_eo);
			spinor_code->copy_to_eoprec_spinorfield_buffer(&in2_eo, sf_in2_eo);
			spinor_code->convert_from_eoprec_device(&in1_eo, &in2_eo, &in1_noneo);

			spinor_code->copy_to_eoprec_spinorfield_buffer(&in3_eo, sf_in3_eo);
			spinor_code->copy_to_eoprec_spinorfield_buffer(&in4_eo, sf_in4_eo);
			spinor_code->convert_from_eoprec_device(&in3_eo, &in4_eo, &in2_noneo);
		}  else {
			in1_noneo.load(sf_in1_noneo);
			in2_noneo.load(sf_in2_noneo);
			spinor_code->convert_to_eoprec_device(&in1_eo, &in2_eo, &in1_noneo);
			spinor_code->convert_to_eoprec_device(&in3_eo, &in4_eo, &in2_noneo);
		}
	}

	hmc_float cpu_back_eo, cpu_back_eo2, cpu_back_eo3, cpu_back_eo4;
	logger.info() << "eo input:";
	logger.info() << "|phi_even_1|^2:";
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in1_eo, &sqnorm);
	sqnorm.dump(&cpu_back_eo);
	logger.info() << cpu_back_eo;
	logger.info() << "|phi_even_2|^2:";
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in2_eo, &sqnorm);
	sqnorm.dump(&cpu_back_eo2);
	logger.info() << cpu_back_eo2;
	logger.info() << "|phi_odd_1|^2:";
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in3_eo, &sqnorm);
	sqnorm.dump(&cpu_back_eo3);
	logger.info() << cpu_back_eo3;
	logger.info() << "|phi_odd_2|^2:";
	spinor_code->set_float_to_global_squarenorm_eoprec_device(&in4_eo, &sqnorm);
	sqnorm.dump(&cpu_back_eo4);
	logger.info() << cpu_back_eo4;

	logger.info() << "run eo force on EVEN and ODD sites...";
	device->fermion_force_eo_device(&in1_eo, &in4_eo, cpu.get_gaugefield(), &out_eo, ODD, params.get_kappa() );
	device->fermion_force_eo_device(&in2_eo, &in3_eo, cpu.get_gaugefield(), &out_eo, EVEN, params.get_kappa() );

	logger.info() << "|force_eo (even) + force_eo (odd)|^2:";
	hmc_float cpu_res_eo;
	gm_code->set_float_to_gaugemomentum_squarenorm_device(&out_eo, &sqnorm);
	sqnorm.dump(&cpu_res_eo);
	logger.info() << cpu_res_eo;

	logger.info() << "non-eo input:";
	hmc_float cpu_back_noneo, cpu_back2_noneo;
	logger.info() << "|phi_1|^2:";
	spinor_code->set_float_to_global_squarenorm_device(&in1_noneo, &sqnorm);
	sqnorm.dump(&cpu_back_noneo);
	logger.info() << cpu_back_noneo;
	logger.info() << "|phi_2|^2:";
	spinor_code->set_float_to_global_squarenorm_device(&in2_noneo, &sqnorm);
	sqnorm.dump(&cpu_back2_noneo);
	logger.info() << cpu_back2_noneo;
	logger.info() << "run noneo force with noneo input...";
	device->fermion_force_device( &in1_noneo, &in2_noneo, cpu.get_gaugefield(), &out_noneo, params.get_kappa());
	logger.info() << "|force_noneo|^2:";
	hmc_float cpu_res_noneo;
	gm_code->set_float_to_gaugemomentum_squarenorm_device(&out_noneo, &sqnorm);
	sqnorm.dump(&cpu_res_noneo);
	logger.info() << cpu_res_noneo;

	logger.info() << "Clear buffers";
	delete[] sf_in1_eo;
	delete[] sf_in2_eo;
	delete[] sf_out_eo;
	delete[] sf_in3_eo;
	delete[] sf_in4_eo;
	delete[] sf_out_noneo;
	delete[] sf_in1_noneo;
	delete[] sf_in2_noneo;

	logger.info() << "Compare eo and non-eo results";
	BOOST_REQUIRE_CLOSE(cpu_res_eo, cpu_res_noneo, params.get_solver_prec() );
	testFloatAgainstInputparameters(cpu_res_eo, params);
	testFloatAgainstInputparameters(cpu_res_noneo, params);
	BOOST_MESSAGE("Test done");
}

BOOST_AUTO_TEST_SUITE( F_FERMION_COMPARE_NONEO_EO )

BOOST_AUTO_TEST_CASE( F_FERMION_COMPARE_NONEO_EO_1 )
{
	test_f_fermion_compare_noneo_eo("/f_fermion_compare_noneo_eo_input_1");
}

BOOST_AUTO_TEST_CASE( F_FERMION_COMPARE_NONEO_EO_2 )
{
	test_f_fermion_compare_noneo_eo("/f_fermion_compare_noneo_eo_input_2");
}

BOOST_AUTO_TEST_CASE( F_FERMION_COMPARE_NONEO_EO_3 )
{
	test_f_fermion_compare_noneo_eo("/f_fermion_compare_noneo_eo_input_3");
}

BOOST_AUTO_TEST_CASE( F_FERMION_COMPARE_NONEO_EO_4 )
{
	test_f_fermion_compare_noneo_eo("/f_fermion_compare_noneo_eo_input_4");
}

BOOST_AUTO_TEST_CASE( F_FERMION_COMPARE_NONEO_EO_5 )
{
	test_f_fermion_compare_noneo_eo("/f_fermion_compare_noneo_eo_input_5");
}

BOOST_AUTO_TEST_CASE( F_FERMION_COMPARE_NONEO_EO_6 )
{
	test_f_fermion_compare_noneo_eo("/f_fermion_compare_noneo_eo_input_6");
}

BOOST_AUTO_TEST_SUITE_END()
