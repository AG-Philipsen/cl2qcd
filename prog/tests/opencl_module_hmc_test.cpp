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

}

void test_gm_squarenorm(std::string inputfile)
{

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
	cl_int err = CL_SUCCESS;
	cl_mem in, sqnorm;
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

	size_t ae_buf_size = device->get_gaugemomentum_buffer_size();
	in = clCreateBuffer(device->get_context(), CL_MEM_READ_WRITE , ae_buf_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	device->importGaugemomentumBuffer(in, reinterpret_cast<ae*>(gm_in));
	sqnorm = clCreateBuffer(device->get_context(), CL_MEM_READ_WRITE, sizeof(hmc_float), 0, &err);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	logger.info() << "|in|^2:";
	device->set_float_to_gaugemomentum_squarenorm_device(in, sqnorm);
	hmc_float cpu_back;
	clEnqueueReadBuffer(device->get_queue(), sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_back, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	logger.info() << cpu_back;

	logger.info() << "Run kernel";
	hmc_float eps = params.get_tau();
	device->md_update_gaugefield_device( in, device->get_gaugefield(), eps);
	logger.info() << "gaugeobservables: ";
	cpu.print_gaugeobservables_from_task(0, 0);

	hmc_float plaq_cpu, tplaq_cpu, splaq_cpu;
	hmc_complex pol_cpu;
	device->gaugeobservables(&plaq_cpu, &tplaq_cpu, &splaq_cpu, &pol_cpu);
	logger.info() << "Finalize device";
	cpu.finalize();

	logger.info() << "Free buffers";
	clReleaseMemObject(in);
	clReleaseMemObject(sqnorm);
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
	cl_int err = CL_SUCCESS;
	cl_mem in, out, sqnorm;
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

	size_t ae_buf_size = device->get_gaugemomentum_buffer_size();
	in = clCreateBuffer(device->get_context(), CL_MEM_READ_WRITE , ae_buf_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	device->importGaugemomentumBuffer(in, reinterpret_cast<ae*>(gm_in));
	out = clCreateBuffer(device->get_context(), CL_MEM_READ_WRITE, ae_buf_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	device->importGaugemomentumBuffer(out, reinterpret_cast<ae*>(gm_out));
	sqnorm = clCreateBuffer(device->get_context(), CL_MEM_READ_WRITE, sizeof(hmc_float), 0, &err);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	logger.info() << "|in|^2:";
	device->set_float_to_gaugemomentum_squarenorm_device(in, sqnorm);
	hmc_float cpu_back;
	clEnqueueReadBuffer(device->get_queue(), sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_back, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	logger.info() << cpu_back;
	logger.info() << "|out|^2:";
	device->set_float_to_gaugemomentum_squarenorm_device(out, sqnorm);
	hmc_float cpu_back2;
	clEnqueueReadBuffer(device->get_queue(), sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_back2, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	logger.info() << cpu_back2;

	logger.info() << "Run kernel";
	hmc_float eps = params.get_tau();
	device->md_update_gaugemomentum_device( in  , out, eps);

	device->set_float_to_gaugemomentum_squarenorm_device(out, sqnorm);
	hmc_float cpu_res;
	clEnqueueReadBuffer(device->get_queue(), sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_res, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	logger.info() << "result:";
	logger.info() << cpu_res;
	logger.info() << "Finalize device";
	cpu.finalize();

	logger.info() << "Free buffers";
	clReleaseMemObject(in);
	clReleaseMemObject(out);
	clReleaseMemObject(sqnorm);
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
	cl_int err = CL_SUCCESS;
	cl_mem out, sqnorm;
	ae * gm_out;

	logger.info() << "create buffers";
	size_t NUM_ELEMENTS_AE = meta::get_vol4d(params) * NDIM * meta::get_su3algebrasize();
	gm_out = new ae[NUM_ELEMENTS_AE];
	fill_with_zero(gm_out, NUM_ELEMENTS_AE);
	BOOST_REQUIRE(gm_out);

	size_t ae_buf_size = device->get_gaugemomentum_buffer_size();
	out = clCreateBuffer(device->get_context(), CL_MEM_READ_WRITE, ae_buf_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	device->importGaugemomentumBuffer(out, reinterpret_cast<ae*>(gm_out));
	sqnorm = clCreateBuffer(device->get_context(), CL_MEM_READ_WRITE, sizeof(hmc_float), 0, &err);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	logger.info() << "|out|^2:";
	device->set_float_to_gaugemomentum_squarenorm_device(out, sqnorm);
	hmc_float cpu_back;
	clEnqueueReadBuffer(device->get_queue(), sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_back, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	logger.info() << cpu_back;

	device->gauge_force_device( device->get_gaugefield(), out);

	logger.info() << "result:";
	hmc_float cpu_res;
	device->set_float_to_gaugemomentum_squarenorm_device(out, sqnorm);
	err = clEnqueueReadBuffer(device->get_queue(), sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_res, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	logger.info() << cpu_res;
	logger.info() << "Finalize device";
	cpu.finalize();

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
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
	cl_int err = CL_SUCCESS;
	cl_mem out, sqnorm;
	ae * gm_out;

	logger.info() << "create buffers";
	size_t NUM_ELEMENTS_AE = meta::get_vol4d(params) * NDIM * meta::get_su3algebrasize();
	gm_out = new ae[NUM_ELEMENTS_AE];
	fill_with_zero(gm_out, NUM_ELEMENTS_AE);
	BOOST_REQUIRE(gm_out);

	size_t ae_buf_size = device->get_gaugemomentum_buffer_size();
	out = clCreateBuffer(device->get_context(), CL_MEM_READ_WRITE, ae_buf_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	device->importGaugemomentumBuffer(out, reinterpret_cast<ae*>(gm_out));
	sqnorm = clCreateBuffer(device->get_context(), CL_MEM_READ_WRITE, sizeof(hmc_float), 0, &err);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	logger.info() << "|out|^2:";
	device->set_float_to_gaugemomentum_squarenorm_device(out, sqnorm);
	hmc_float cpu_back;
	clEnqueueReadBuffer(device->get_queue(), sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_back, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	logger.info() << cpu_back;

	device->gauge_force_tlsym_device( device->get_gaugefield(), out);

	logger.info() << "result:";
	hmc_float cpu_res;
	device->set_float_to_gaugemomentum_squarenorm_device(out, sqnorm);
	err = clEnqueueReadBuffer(device->get_queue(), sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_res, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	logger.info() << cpu_res;
	logger.info() << "Finalize device";
	cpu.finalize();

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
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
	cl_int err = CL_SUCCESS;
	cl_mem out, sqnorm;
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

	size_t sf_buf_size = meta::get_spinorfieldsize(params) * sizeof(spinor);

	const Plain<spinor> in1(NUM_ELEMENTS_SF, device->get_device());
	const Plain<spinor> in2(NUM_ELEMENTS_SF, device->get_device());
	in1.load(sf_in1);
	in2.load(sf_in2);
	out = clCreateBuffer(device->get_context(), CL_MEM_READ_WRITE, device->get_gaugemomentum_buffer_size(), 0, &err);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	device->importGaugemomentumBuffer(out, ae_out);
	sqnorm = clCreateBuffer(device->get_context(), CL_MEM_READ_WRITE, sizeof(hmc_float), 0, &err);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	logger.info() << "|phi_1|^2:";
	hmc_float cpu_back;
	device->set_float_to_global_squarenorm_device(&in1, sqnorm);
	err = clEnqueueReadBuffer(device->get_queue(), sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_back, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	logger.info() << cpu_back;
	logger.info() << "|phi_2|^2:";
	hmc_float cpu_back2;
	device->set_float_to_global_squarenorm_device(&in2, sqnorm);
	err = clEnqueueReadBuffer(device->get_queue(), sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_back2, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	logger.info() << cpu_back2;
	logger.info() << "Run kernel";
	device->fermion_force_device( &in1, &in2, device->get_gaugefield()  , out, params.get_kappa());
	logger.info() << "result:";
	hmc_float cpu_res;
	device->set_float_to_gaugemomentum_squarenorm_device(out, sqnorm);
	err = clEnqueueReadBuffer(device->get_queue(), sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_res, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	logger.info() << cpu_res;
	logger.info() << "Finalize device";
	cpu.finalize();

	logger.info() << "Clear buffers";
	clReleaseMemObject(out);
	clReleaseMemObject(sqnorm);
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
	cl_int err = CL_SUCCESS;
	cl_mem out;
	cl_mem sqnorm;
	spinor * sf_in1;
	spinor * sf_in2;
	spinor * sf_in3;
	spinor * sf_in4;
	ae * ae_out;

	logger.info() << "fill buffers";
	size_t NUM_ELEMENTS_SF =  meta::get_eoprec_spinorfieldsize(params);
	size_t NUM_ELEMENTS_AE = meta::get_vol4d(params) * NDIM * meta::get_su3algebrasize();

	sf_in1 = new spinor[NUM_ELEMENTS_SF];
	sf_in2 = new spinor[NUM_ELEMENTS_SF];
	sf_in3 = new spinor[NUM_ELEMENTS_SF];
	sf_in4 = new spinor[NUM_ELEMENTS_SF];
	ae_out = new ae[NUM_ELEMENTS_AE];

	//use the variable use_cg to switch between cold and random input sf
	if(params.get_solver() == meta::Inputparameters::cg) {
		fill_sf_with_one(sf_in1, NUM_ELEMENTS_SF);
		fill_sf_with_one(sf_in2, NUM_ELEMENTS_SF);
		fill_sf_with_one(sf_in3, NUM_ELEMENTS_SF);
		fill_sf_with_one(sf_in4, NUM_ELEMENTS_SF);
	} else {
		fill_sf_with_random_eo(sf_in1, sf_in2, NUM_ELEMENTS_SF, 123456);
		fill_sf_with_random_eo(sf_in3, sf_in4, NUM_ELEMENTS_SF, 789101);
	}
	fill_with_zero(ae_out, NUM_ELEMENTS_AE);
	BOOST_REQUIRE(sf_in1);
	BOOST_REQUIRE(sf_in2);
	BOOST_REQUIRE(sf_in3);
	BOOST_REQUIRE(sf_in4);
	BOOST_REQUIRE(ae_out);

	const Spinor in1(NUM_ELEMENTS_SF, device->get_device());
	const Spinor in2(NUM_ELEMENTS_SF, device->get_device());
	const Spinor in3(NUM_ELEMENTS_SF, device->get_device());
	const Spinor in4(NUM_ELEMENTS_SF, device->get_device());
	out = clCreateBuffer(device->get_context(), CL_MEM_WRITE_ONLY, device->get_gaugemomentum_buffer_size(), 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	sqnorm = clCreateBuffer(device->get_context(), CL_MEM_READ_WRITE, sizeof(hmc_float), 0, &err);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	device->copy_to_eoprec_spinorfield_buffer(&in1, sf_in1);
	device->copy_to_eoprec_spinorfield_buffer(&in2, sf_in2);
	device->copy_to_eoprec_spinorfield_buffer(&in3, sf_in3);
	device->copy_to_eoprec_spinorfield_buffer(&in4, sf_in4);
	device->importGaugemomentumBuffer(out, ae_out);

	//switch according to "use_pointsource"
	hmc_float cpu_res, cpu_back, cpu_back2, cpu_back3, cpu_back4;
	//interprete Y = (in1, in2) X = (in3, in4)
	//Y_odd = in2, Y_even = in1, X_odd = in4, X_even = in3

	if(params.get_use_pointsource()) {
		logger.info() << "|phi_even_1|^2:";
		device->set_float_to_global_squarenorm_eoprec_device(&in1, sqnorm);
		err = clEnqueueReadBuffer(device->get_queue(), sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_back, 0, 0, 0);
		BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
		logger.info() << cpu_back;
		logger.info() << "|phi_even_2|^2:";
		device->set_float_to_global_squarenorm_eoprec_device(&in4, sqnorm);
		err = clEnqueueReadBuffer(device->get_queue(), sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_back2, 0, 0, 0);
		BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
		logger.info() << cpu_back2;
		logger.info() << "Run kernel";
		int tmp = EVEN;
		//this is then force(Y_even, X_odd) == force(in1, in4)
		device->fermion_force_eo_device(&in1, &in4, device->get_gaugefield(), out, tmp, params.get_kappa() );

		logger.info() << "|force (even)|^2:";
	} else {
		logger.info() << "|phi_odd_1|^2:";
		device->set_float_to_global_squarenorm_eoprec_device(&in2, sqnorm);
		err = clEnqueueReadBuffer(device->get_queue(), sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_back3, 0, 0, 0);
		BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
		logger.info() << cpu_back3;
		logger.info() << "|phi_odd_2|^2:";
		device->set_float_to_global_squarenorm_eoprec_device(&in3, sqnorm);
		err = clEnqueueReadBuffer(device->get_queue(), sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_back4, 0, 0, 0);
		BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
		logger.info() << cpu_back4;
		logger.info() << "Run kernel";

		int tmp = ODD;
		//this is then force(Y_odd, X_even) == force(in2, in3)
		device->fermion_force_eo_device(&in2, &in3, device->get_gaugefield(), out, tmp, params.get_kappa() );
		logger.info() << "|force (odd)|^2:";
	}

	device->set_float_to_gaugemomentum_squarenorm_device(out, sqnorm);
	err = clEnqueueReadBuffer(device->get_queue(), sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_res, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	logger.info() << cpu_res;
	logger.info() << "Finalize device";
	cpu.finalize();

	logger.info() << "Clear buffers";
	clReleaseMemObject(out);
	clReleaseMemObject(sqnorm);

	delete[] sf_in1;
	delete[] sf_in2;
	delete[] sf_in3;
	delete[] sf_in4;
	delete[] ae_out;

	testFloatAgainstInputparameters(cpu_res, params);
	BOOST_MESSAGE("Test done");
}

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

BOOST_AUTO_TEST_SUITE( SET_ZERO_GAUGEMOMENTUM )

BOOST_AUTO_TEST_CASE( SET_ZERO_GAUGEMOMENTUM_1 )
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!!");
	test_set_zero_gm("/set_zero_gaugemomentum_input_1");
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
	BOOST_MESSAGE("NOT YET IMPLEMENTED!!");
	test_gm_squarenorm("/gm_squarenorm_input_1");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( GF_UPDATE )

BOOST_AUTO_TEST_CASE( GF_UPDATE_1 )
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!!");
	test_gf_update("/gf_update_input_1");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_UPDATE )

BOOST_AUTO_TEST_CASE( F_UPDATE_1 )
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!!");
	test_f_update("/f_update_input_1");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_GAUGE )

BOOST_AUTO_TEST_CASE( F_GAUGE_1 )
{
	test_f_gauge("/f_gauge_input_1");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_GAUGE_TLSYM )

BOOST_AUTO_TEST_CASE( F_GAUGE_TLSYM_1 )
{
	test_f_gauge_tlsym("/f_gauge_tlsym_input_1");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_FERMION )

BOOST_AUTO_TEST_CASE( F_FERMION_1 )
{
	test_f_fermion("/f_fermion_input_1");
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

BOOST_AUTO_TEST_SUITE_END()
