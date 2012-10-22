#include "../opencl_module_hmc.h"
#include "../gaugefield_hybrid.h"

#include "../meta/util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE_HMC
#include <boost/test/unit_test.hpp>

#include "test_util.h"


extern std::string const version;
std::string const version = "0.1";

class TestGaugefield : public Gaugefield_hybrid {

public:
	TestGaugefield(const hardware::System * system) : Gaugefield_hybrid(system) {
		auto inputfile = system->get_inputparameters();
		init(1, inputfile.get_use_gpu() ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU);
		meta::print_info_hmc("test program", inputfile);
		print_gaugeobservables(0);
	};

	virtual void init_tasks();
	virtual void finalize_opencl();

	Opencl_Module_Hmc * get_device();
};

void TestGaugefield::init_tasks()
{
	opencl_modules = new Opencl_Module* [get_num_tasks()];
	meta::Counter counter1, counter2, counter3, counter4;
	opencl_modules[0] = new Opencl_Module_Hmc(get_parameters(), get_device_for_task(0), &counter1, &counter2, &counter3, &counter4);
	opencl_modules[0]->init();
}

void TestGaugefield::finalize_opencl()
{
	Gaugefield_hybrid::finalize_opencl();
}

Opencl_Module_Hmc* TestGaugefield::get_device()
{
	return static_cast<Opencl_Module_Hmc*>(opencl_modules[0]);
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
	for(int i = 0; i < size; ++i) {
		prng_init(seed);
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
	//  Random rnd_loc(seed);
	for(int i = 0; i < size; ++i) {
		prng_init(seed);
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
  logger.info() << "build opencl_module_hmc";
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	logger.info() << "Finalize device";
	cpu.finalize();
	BOOST_MESSAGE("Test done");
}

void test_generate_gaussian_spinorfield(std::string inputfile)
{

}

void test_generate_gaussian_spinorfield_eo(std::string inputfile)
{

}

void test_generate_gaussian_gaugemomenta(std::string inputfile)
{

}

void test_stout_smear_fermion_force(std::string inputfile)
{

}

void test_set_zero_gm(std::string inputfile)
{
	std::string kernelName = "set_zero_gm";
	printKernelInfo(kernelName);

	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	Opencl_Module_Hmc * device = cpu.get_device();
	hmc_float * gm_in;

	logger.info() << "create buffers";
	size_t NUM_ELEMENTS_AE = meta::get_vol4d(params) * NDIM * meta::get_su3algebrasize();
	gm_in = new hmc_float[NUM_ELEMENTS_AE];
	fill_with_random(gm_in, NUM_ELEMENTS_AE, 123456);
	BOOST_REQUIRE(gm_in);

	hardware::buffers::Gaugemomentum in(meta::get_vol4d(params) * NDIM, device->get_device());
	device->importGaugemomentumBuffer(&in, reinterpret_cast<ae*>(gm_in));
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	logger.info() << "|in|^2:";
	device->set_float_to_gaugemomentum_squarenorm_device(&in, &sqnorm);
	hmc_float cpu_back;
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;

	logger.info() << "Run kernel";
	device->set_zero_gaugemomentum(&in);

	logger.info() << "|out|^2:";
	device->set_float_to_gaugemomentum_squarenorm_device(&in, &sqnorm);
	hmc_float cpu_back2;
	sqnorm.dump(&cpu_back2);
	logger.info() << cpu_back2;

	logger.info() << "Finalize device";
	cpu.finalize();

	logger.info() << "Free buffers";
	delete[] gm_in;

	testFloatAgainstInputparameters(cpu_back2, params);
	BOOST_MESSAGE("Test done");

}

void test_gm_squarenorm(std::string inputfile)
{
	std::string kernelName = "gaugemomenta squarenorm";
	printKernelInfo(kernelName);

	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	Opencl_Module_Hmc * device = cpu.get_device();
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

	hardware::buffers::Gaugemomentum in(meta::get_vol4d(params) * NDIM, device->get_device());
	device->importGaugemomentumBuffer(&in, reinterpret_cast<ae*>(gm_in));
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	logger.info() << "Run kernel";
	logger.info() << "|in|^2:";
	device->set_float_to_gaugemomentum_squarenorm_device(&in, &sqnorm);
	hmc_float cpu_back;
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;

	logger.info() << "Finalize device";
	cpu.finalize();

	logger.info() << "Free buffers";
	delete[] gm_in;

	testFloatAgainstInputparameters(cpu_back, params);
	BOOST_MESSAGE("Test done");
}

void test_gm_convert_to_soa(std::string inputfile)
{

}

void test_gm_convert_from_soa(std::string inputfile)
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
	Opencl_Module_Hmc * device = cpu.get_device();
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

	hardware::buffers::Gaugemomentum in(meta::get_vol4d(params) * NDIM, device->get_device());
	device->importGaugemomentumBuffer(&in, reinterpret_cast<ae*>(gm_in));
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	logger.info() << "|in|^2:";
	device->set_float_to_gaugemomentum_squarenorm_device(&in, &sqnorm);
	hmc_float cpu_back;
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;

	logger.info() << "Run kernel";
	hmc_float eps = params.get_tau();
	device->md_update_gaugefield_device(&in, device->get_gaugefield(), eps);
	logger.info() << "gaugeobservables: ";
	cpu.print_gaugeobservables_from_task(0, 0);

	hmc_float plaq_cpu, tplaq_cpu, splaq_cpu;
	hmc_complex pol_cpu;
	device->gaugeobservables(&plaq_cpu, &tplaq_cpu, &splaq_cpu, &pol_cpu);
	logger.info() << "Finalize device";
	cpu.finalize();

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
	Opencl_Module_Hmc * device = cpu.get_device();
	hmc_float * gm_in;
	hmc_float * gm_out;

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
	device->importGaugemomentumBuffer(&in, reinterpret_cast<ae*>(gm_in));
	hardware::buffers::Gaugemomentum out(meta::get_vol4d(params) * NDIM, device->get_device());
	device->importGaugemomentumBuffer(&out, reinterpret_cast<ae*>(gm_out));
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	logger.info() << "|in|^2:";
	device->set_float_to_gaugemomentum_squarenorm_device(&in, &sqnorm);
	hmc_float cpu_back;
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;
	logger.info() << "|out|^2:";
	device->set_float_to_gaugemomentum_squarenorm_device(&out, &sqnorm);
	hmc_float cpu_back2;
	sqnorm.dump(&cpu_back2);
	logger.info() << cpu_back2;

	logger.info() << "Run kernel";
	hmc_float eps = params.get_tau();
	device->md_update_gaugemomentum_device(&in, &out, eps);

	device->set_float_to_gaugemomentum_squarenorm_device(&out, &sqnorm);
	hmc_float cpu_res;
	sqnorm.dump(&cpu_res);
	logger.info() << "result:";
	logger.info() << cpu_res;
	logger.info() << "Finalize device";
	cpu.finalize();

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
	Opencl_Module_Hmc * device = cpu.get_device();
	ae * gm_out;

	logger.info() << "create buffers";
	size_t NUM_ELEMENTS_AE = meta::get_vol4d(params) * NDIM * meta::get_su3algebrasize();
	gm_out = new ae[NUM_ELEMENTS_AE];
	fill_with_zero(gm_out, NUM_ELEMENTS_AE);
	BOOST_REQUIRE(gm_out);

	hardware::buffers::Gaugemomentum out(meta::get_vol4d(params) * NDIM, device->get_device());
	device->importGaugemomentumBuffer(&out, reinterpret_cast<ae*>(gm_out));
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	logger.info() << "|out|^2:";
	device->set_float_to_gaugemomentum_squarenorm_device(&out, &sqnorm);
	hmc_float cpu_back;
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;

	device->gauge_force_device( device->get_gaugefield(), &out);

	logger.info() << "result:";
	hmc_float cpu_res;
	device->set_float_to_gaugemomentum_squarenorm_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");

	logger.info() << "Finalize device";
	delete[] gm_out;
	cpu.finalize();
}

void test_f_gauge_tlsym(std::string inputfile)
{
	std::string kernelName = "f_gauge_tlsym";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	Opencl_Module_Hmc * device = cpu.get_device();
	ae * gm_out;

	logger.info() << "create buffers";
	size_t NUM_ELEMENTS_AE = meta::get_vol4d(params) * NDIM * meta::get_su3algebrasize();
	gm_out = new ae[NUM_ELEMENTS_AE];
	fill_with_zero(gm_out, NUM_ELEMENTS_AE);
	BOOST_REQUIRE(gm_out);

	hardware::buffers::Gaugemomentum out(meta::get_vol4d(params) * NDIM, device->get_device());
	device->importGaugemomentumBuffer(&out, reinterpret_cast<ae*>(gm_out));
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	logger.info() << "|out|^2:";
	device->set_float_to_gaugemomentum_squarenorm_device(&out, &sqnorm);
	hmc_float cpu_back;
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;

	device->gauge_force_tlsym_device( device->get_gaugefield(), &out);

	logger.info() << "result:";
	hmc_float cpu_res;
	device->set_float_to_gaugemomentum_squarenorm_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;
	logger.info() << "Finalize device";
	cpu.finalize();

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
	Opencl_Module_Hmc * device = cpu.get_device();
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

	const Plain<spinor> in1(NUM_ELEMENTS_SF, device->get_device());
	const Plain<spinor> in2(NUM_ELEMENTS_SF, device->get_device());
	in1.load(sf_in1);
	in2.load(sf_in2);
	Gaugemomentum out(meta::get_vol4d(params) * NDIM, device->get_device());
	device->importGaugemomentumBuffer(&out, reinterpret_cast<ae*>(ae_out));
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());

	logger.info() << "|phi_1|^2:";
	hmc_float cpu_back;
	device->set_float_to_global_squarenorm_device(&in1, &sqnorm);
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;
	logger.info() << "|phi_2|^2:";
	hmc_float cpu_back2;
	device->set_float_to_global_squarenorm_device(&in2, &sqnorm);
	sqnorm.dump(&cpu_back2);
	logger.info() << cpu_back2;
	logger.info() << "Run kernel";
	device->fermion_force_device( &in1, &in2, device->get_gaugefield(), &out, params.get_kappa());
	logger.info() << "result:";
	hmc_float cpu_res;
	device->set_float_to_gaugemomentum_squarenorm_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;
	logger.info() << "Finalize device";
	cpu.finalize();

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
	Opencl_Module_Hmc * device = cpu.get_device();
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

	const Spinor in1(NUM_ELEMENTS_SF, device->get_device());
	const Spinor in2(NUM_ELEMENTS_SF, device->get_device());
	Gaugemomentum out(meta::get_vol4d(params) * NDIM, device->get_device());
	hardware::buffers::Plain<hmc_float> sqnorm(1, device->get_device());
	device->copy_to_eoprec_spinorfield_buffer(&in1, sf_in1);
	device->copy_to_eoprec_spinorfield_buffer(&in2, sf_in2);
	device->importGaugemomentumBuffer(&out, ae_out);


	hmc_float cpu_res, cpu_back, cpu_back2;
	logger.info() << "|in_1|^2:";
	device->set_float_to_global_squarenorm_eoprec_device(&in1, &sqnorm);
	sqnorm.dump(&cpu_back);
	logger.info() << cpu_back;
	logger.info() << "|in_2|^2:";
	device->set_float_to_global_squarenorm_eoprec_device(&in2, &sqnorm);
	sqnorm.dump(&cpu_back2);
	logger.info() << cpu_back2;
	logger.info() << "Run kernel";

	//switch according to "use_pointsource"
	if(params.get_use_pointsource()) {
		int tmp = EVEN;
		device->fermion_force_eo_device(&in1, &in2, device->get_gaugefield(), &out, tmp, params.get_kappa() );
	} else {
		int tmp = ODD;
		device->fermion_force_eo_device(&in1, &in2, device->get_gaugefield(), &out, tmp, params.get_kappa() );
	}
	logger.info() << "|force|^2:";
	device->set_float_to_gaugemomentum_squarenorm_device(&out, &sqnorm);
	sqnorm.dump(&cpu_res);
	logger.info() << cpu_res;
	logger.info() << "Finalize device";
	cpu.finalize();

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
  test_build("/opencl_module_hmc_build_input_1");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(GENERATE_GAUSSIAN_SPINORFIELD  )

BOOST_AUTO_TEST_CASE( GENERATE_GAUSSIAN_SPINORFIELD_1 )
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!!");
	test_generate_gaussian_spinorfield("/generate_gaussian_spinorfield_input_1");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(GENERATE_GAUSSIAN_SPINORFIELD_EO  )

BOOST_AUTO_TEST_CASE( GENERATE_GAUSSIAN_SPINORFIELD_EO_1 )
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!!");
	test_generate_gaussian_spinorfield_eo("/generate_gaussian_spinorfield_eo_input_1");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(GENERATE_GAUSSIAN_GAUGEMOMENTA  )

BOOST_AUTO_TEST_CASE(GENERATE_GAUSSIAN_GAUGEMOMENTA_1 )
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!!");
	test_generate_gaussian_gaugemomenta("/generate_gaussian_gaugemomenta_input_1");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( STOUT_SMEAR_FERMION_FORCE )

BOOST_AUTO_TEST_CASE(STOUT_SMEAR_FERMION_FORCE_1 )
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!!");
	test_stout_smear_fermion_force("/stout_smear_fermion_force_input_1");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( GM_CONVERT_TO_SOA )

BOOST_AUTO_TEST_CASE( GM_CONVERT_TO_SOA_1 )
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!!");
	test_gm_convert_to_soa("/gm_convert_to_soa_input_1");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( GM_CONVERT_FROM_SOA )

BOOST_AUTO_TEST_CASE( GM_CONVERT_FROM_SOA_1 )
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!!");
	test_gm_convert_from_soa("/gm_convert_from_soa_input_1");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( GM_SQUARENORM )

BOOST_AUTO_TEST_CASE(GM_SQUARENORM_1  )
{
	test_gm_squarenorm("/gm_squarenorm_input_1");
}

BOOST_AUTO_TEST_CASE(GM_SQUARENORM_2  )
{
	test_gm_squarenorm("/gm_squarenorm_input_2");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( SET_ZERO_GM )

BOOST_AUTO_TEST_CASE( SET_ZERO_GM_1 )
{
	test_set_zero_gm("/gm_set_zero_input_1");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( GF_UPDATE )

BOOST_AUTO_TEST_CASE(GF_UPDATE_1 )
{
	test_gf_update("/gf_update_input_1");
}

BOOST_AUTO_TEST_CASE( GF_UPDATE_2 )
{
  BOOST_MESSAGE("THIS TEST HAS TO BE INVESTIGATED!! IT COULD HINT TO AN ERROR IN THE FUNCTION!");
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

BOOST_AUTO_TEST_SUITE_END()



