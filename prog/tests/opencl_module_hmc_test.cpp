#include "../opencl_module_hmc.h"
#include "../gaugefield_hybrid.h"

#include "../meta/util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE F_fermion_test
#include <boost/test/unit_test.hpp>

#include "test_util.h"


extern std::string const version;
std::string const version = "0.1";
std::string const exec_name = "f_fermion";

class TestGaugefield : public Gaugefield_hybrid {

public:
	TestGaugefield(const hardware::System * system) : Gaugefield_hybrid(system) {
		auto inputfile = system->get_inputparameters();
		init(1, inputfile.get_use_gpu() ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU);
		meta::print_info_hmc("test program", inputfile);
	};

	virtual void init_tasks();
	virtual void finalize_opencl();

  cl_mem in1, in2, in3, in4, out;
  cl_mem in1_eo, in2_eo, in3_eo, in4_eo, out_eo;
	cl_mem sqnorm;
  Opencl_Module_Hmc * get_device();
private:
	void fill_buffers();
	void clear_buffers();
	spinor * sf_in;
  spinor * sf_in1;
  spinor * sf_in2;
  spinor * sf_in1_eo;
  spinor * sf_in2_eo;
  spinor * sf_in3_eo;
  spinor * sf_in4_eo;
  spinor * sf_out;

	ae * ae_out;
};

void TestGaugefield::init_tasks()
{
	opencl_modules = new Opencl_Module* [get_num_tasks()];
	meta::Counter counter1, counter2, counter3, counter4;
	opencl_modules[0] = new Opencl_Module_Hmc(get_parameters(), get_device_for_task(0), &counter1, &counter2, &counter3, &counter4);
	opencl_modules[0]->init();

	fill_buffers();
}

void TestGaugefield::finalize_opencl()
{
	clear_buffers();
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

void fill_with_one(hmc_float * sf_in, int size)
{
	for(int i = 0; i < size; ++i) {
		sf_in[i] = 1.;
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

void TestGaugefield::fill_buffers()
{
	// don't invoke parent function as we don't require the original buffers
	cl_int err;
	cl_context context = opencl_modules[0]->get_context();
	int NUM_ELEMENTS_SF = meta::get_spinorfieldsize(get_parameters());
	int NUM_ELEMENTS_SF_EO =  meta::get_eoprec_spinorfieldsize(get_parameters());
	int NUM_ELEMENTS_AE = meta::get_vol4d(get_parameters()) * NDIM * meta::get_su3algebrasize();

	sf_in1 = new spinor[NUM_ELEMENTS_SF];
	sf_in2 = new spinor[NUM_ELEMENTS_SF];

	sf_in1_eo = new spinor[NUM_ELEMENTS_SF];
	sf_in2_eo = new spinor[NUM_ELEMENTS_SF];
	sf_in3_eo = new spinor[NUM_ELEMENTS_SF];
	sf_in4_eo = new spinor[NUM_ELEMENTS_SF];

	ae_out = new ae[NUM_ELEMENTS_AE];

	//use the variable use_cg to switch between cold and random input sf
	if(get_parameters().get_solver() == meta::Inputparameters::cg) {
		fill_sf_with_one(sf_in1, NUM_ELEMENTS_SF);
		fill_sf_with_one(sf_in2, NUM_ELEMENTS_SF);
	} else {
		fill_sf_with_random(sf_in1, NUM_ELEMENTS_SF, 123456);
		fill_sf_with_random(sf_in2, NUM_ELEMENTS_SF, 789101);
	}
	BOOST_REQUIRE(sf_in1);
	BOOST_REQUIRE(sf_in2);

	//use the variable use_cg to switch between cold and random input sf
	if(get_parameters().get_solver() == meta::Inputparameters::cg) {
		fill_sf_with_one(sf_in1_eo, NUM_ELEMENTS_SF_EO);
		fill_sf_with_one(sf_in2_eo, NUM_ELEMENTS_SF_EO);
		fill_sf_with_one(sf_in3_eo, NUM_ELEMENTS_SF_EO);
		fill_sf_with_one(sf_in4_eo, NUM_ELEMENTS_SF_EO);
	} else {
		fill_sf_with_random_eo(sf_in1_eo, sf_in2_eo, NUM_ELEMENTS_SF_EO, 123456);
		fill_sf_with_random_eo(sf_in3_eo, sf_in4_eo, NUM_ELEMENTS_SF_EO, 789101);
	}
	BOOST_REQUIRE(sf_in1_eo);
	BOOST_REQUIRE(sf_in2_eo);
	BOOST_REQUIRE(sf_in3_eo);
	BOOST_REQUIRE(sf_in4_eo);

	fill_with_zero(ae_out, NUM_ELEMENTS_AE);

	size_t sf_buf_size = meta::get_spinorfieldsize(get_parameters()) * sizeof(spinor);
	//create buffer for sf on device (and copy sf_in to both for convenience)

	Opencl_Module_Hmc * device = this->get_device();
	size_t sf_eoprec_buffer_size = device->get_eoprec_spinorfield_buffer_size();

	in1 = clCreateBuffer(context, CL_MEM_READ_ONLY , sf_buf_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	in2 = clCreateBuffer(context, CL_MEM_READ_ONLY , sf_buf_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	err = clEnqueueWriteBuffer(device->get_queue(), in1, CL_TRUE, 0, sf_buf_size, sf_in1, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	err = clEnqueueWriteBuffer(device->get_queue(), in2, CL_TRUE, 0, sf_buf_size, sf_in2, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	out = clCreateBuffer(context, CL_MEM_READ_WRITE, device->get_gaugemomentum_buffer_size(), 0, &err);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	device->importGaugemomentumBuffer(out, reinterpret_cast<ae*>(ae_out));
	/*
	if(get_parameters().get_use_eo() ) {
	  size_t eo_buf_size = device->get_eoprec_spinorfield_buffer_size();
	  in1_eo = clCreateBuffer(context, CL_MEM_READ_WRITE, eo_buf_size, 0, &err );
	  BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	  in2_eo = clCreateBuffer(context, CL_MEM_READ_WRITE, eo_buf_size, 0, &err );
	  BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	  out_eo = clCreateBuffer(context, CL_MEM_READ_WRITE, eo_buf_size, 0, &err );
	  BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	  device->convert_to_eoprec_device(in1_eo, in2_eo, in1);
	} else {
	  in1_eo = 0;
	  in2_eo = 0;
	  out_eo = 0;
	}
	*/
	in1_eo = clCreateBuffer(context, CL_MEM_READ_ONLY , sf_eoprec_buffer_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	in2_eo = clCreateBuffer(context, CL_MEM_READ_ONLY , sf_eoprec_buffer_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	in3_eo = clCreateBuffer(context, CL_MEM_READ_ONLY , sf_eoprec_buffer_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	in4_eo = clCreateBuffer(context, CL_MEM_READ_ONLY , sf_eoprec_buffer_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	err = clEnqueueWriteBuffer(device->get_queue(), in1_eo, CL_TRUE, 0, sf_eoprec_buffer_size, sf_in1_eo, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	err = clEnqueueWriteBuffer(device->get_queue(), in2_eo, CL_TRUE, 0, sf_eoprec_buffer_size, sf_in2_eo, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	err = clEnqueueWriteBuffer(device->get_queue(), in3_eo, CL_TRUE, 0, sf_eoprec_buffer_size, sf_in3_eo, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	err = clEnqueueWriteBuffer(device->get_queue(), in4_eo, CL_TRUE, 0, sf_eoprec_buffer_size, sf_in4_eo, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	//create buffer for squarenorm on device
	sqnorm = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float), 0, &err);
}

void TestGaugefield::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers
  	//clReleaseMemObject(in1);
	//clReleaseMemObject(in2);
	//clReleaseMemObject(out);
	//clReleaseMemObject(in1_eo);
	//clReleaseMemObject(in2_eo);
	//clReleaseMemObject(in3_eo);
	//clReleaseMemObject(in4_eo);
	//clReleaseMemObject(out_eo);
	//clReleaseMemObject(sqnorm);
  /*
	delete[] sf_in1;
	delete[] sf_in2;
	delete[] sf_in1_eo;
	delete[] sf_in2_eo;
	delete[] sf_in3_eo;
	delete[] sf_in4_eo;
	delete[] sf_out;
	delete[] ae_out;
  */
}

Opencl_Module_Hmc* TestGaugefield::get_device()
{
	return static_cast<Opencl_Module_Hmc*>(opencl_modules[0]);
}

void test_f_gauge(std::string inputfile)
{
  std::string kernelName = "f_fermion";
  printKernelInfo(kernelName);
  logger.info() << "Init device";
  meta::Inputparameters params = create_parameters(inputfile);
  hardware::System system(params);
  TestGaugefield cpu(&system);
  Opencl_Module_Hmc * device = cpu.get_device();
  cl_int err = CL_SUCCESS;

  device->gauge_force_device( device->get_gaugefield(), cpu.out);

  logger.info() << "result:";
  hmc_float cpu_res;
  device->set_float_to_gaugemomentum_squarenorm_device(cpu.out, cpu.sqnorm);
  err = clEnqueueReadBuffer(device->get_queue(), cpu.sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_res, 0, 0, 0);
  BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
  logger.info() << cpu_res;
  logger.info() << "Finalize device";
  cpu.finalize();

  logger.info() << "Choosing reference value and acceptance precision";
  hmc_float ref_val = params.get_test_ref_value();
  logger.info() << "reference value:\t" << ref_val;
  hmc_float prec = params.get_solver_prec();
  logger.info() << "acceptance precision: " << prec;
  
  logger.info() << "Compare result to reference value";
  BOOST_REQUIRE_CLOSE(cpu_res, ref_val, prec);
  BOOST_MESSAGE("Test done");
}

void test_f_gauge_tlsym(std::string inputfile)
{
  std::string kernelName = "f_fermion";
  printKernelInfo(kernelName);
  logger.info() << "Init device";
  meta::Inputparameters params = create_parameters(inputfile);
  hardware::System system(params);
  TestGaugefield cpu(&system);
  Opencl_Module_Hmc * device = cpu.get_device();
  cl_int err = CL_SUCCESS;

  device->gauge_force_tlsym_device( device->get_gaugefield(), cpu.out);

  logger.info() << "result:";
  hmc_float cpu_res;
  device->set_float_to_gaugemomentum_squarenorm_device(cpu.out, cpu.sqnorm);
  err = clEnqueueReadBuffer(device->get_queue(), cpu.sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_res, 0, 0, 0);
  BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
  logger.info() << cpu_res;
  logger.info() << "Finalize device";
  cpu.finalize();

  logger.info() << "Choosing reference value and acceptance precision";
  hmc_float ref_val = params.get_test_ref_value();
  logger.info() << "reference value:\t" << ref_val;
  hmc_float prec = params.get_solver_prec();
  logger.info() << "acceptance precision: " << prec;
  
  logger.info() << "Compare result to reference value";
  BOOST_REQUIRE_CLOSE(cpu_res, ref_val, prec);
  BOOST_MESSAGE("Test done");
}

void test_f_fermion(std::string inputfile)
{
  std::string kernelName = "f_fermion";
  printKernelInfo(kernelName);
  logger.info() << "Init device";
  meta::Inputparameters params = create_parameters(inputfile);
  hardware::System system(params);
  TestGaugefield cpu(&system);
  Opencl_Module_Hmc * device = cpu.get_device();
  cl_int err = CL_SUCCESS;
  logger.info() << "|phi_1|^2:";
  hmc_float cpu_back;
  device->set_float_to_global_squarenorm_device(cpu.in1, cpu.sqnorm);
  err = clEnqueueReadBuffer(device->get_queue(), cpu.sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_back, 0, 0, 0);
  BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
  logger.info() << cpu_back;
  logger.info() << "|phi_2|^2:";
  hmc_float cpu_back2;
  device->set_float_to_global_squarenorm_device(cpu.in2, cpu.sqnorm);
  err = clEnqueueReadBuffer(device->get_queue(), cpu.sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_back2, 0, 0, 0);
  BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
  logger.info() << cpu_back2;
  logger.info() << "Run kernel";
  device->fermion_force_device( cpu.in1, cpu.in2, device->get_gaugefield()  , cpu.out, params.get_kappa());
  logger.info() << "result:";
  hmc_float cpu_res;
  device->set_float_to_gaugemomentum_squarenorm_device(cpu.out, cpu.sqnorm);
  err = clEnqueueReadBuffer(device->get_queue(), cpu.sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_res, 0, 0, 0);
  BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
  logger.info() << cpu_res;
  logger.info() << "Finalize device";
  cpu.finalize();

  logger.info() << "Choosing reference value and acceptance precision";
  hmc_float ref_val = params.get_test_ref_value();
  logger.info() << "reference value:\t" << ref_val;
  hmc_float prec = params.get_solver_prec();
  logger.info() << "acceptance precision: " << prec;
  
  logger.info() << "Compare result to reference value";
  BOOST_REQUIRE_CLOSE(cpu_res, ref_val, prec);
  BOOST_MESSAGE("Test done");
}

void test_f_fermion_eo(std::string inputfile)
{
  std::string kernelName = "f_fermion_eo";
  printKernelInfo(kernelName);
  logger.info() << "Init device";
  meta::Inputparameters params = create_parameters(inputfile);
  hardware::System system(params);
  TestGaugefield cpu(&system);
  Opencl_Module_Hmc * device = cpu.get_device();
  cl_int err = CL_SUCCESS;

  //switch according to "use_pointsource"
  hmc_float cpu_res, cpu_back, cpu_back2, cpu_back3, cpu_back4;
  //interprete Y = (in1, in2) X = (in3, in4)
  //Y_odd = in2, Y_even = in1, X_odd = in4, X_even = in3
  
  if(params.get_use_pointsource()) {
    logger.info() << "|phi_even_1|^2:";
    device->set_float_to_global_squarenorm_eoprec_device(cpu.in1_eo, cpu.sqnorm);
    err = clEnqueueReadBuffer(device->get_queue(), cpu.sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_back, 0, 0, 0);
    BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
    logger.info() << cpu_back;
    logger.info() << "|phi_even_2|^2:";
    device->set_float_to_global_squarenorm_eoprec_device(cpu.in2_eo, cpu.sqnorm);
    err = clEnqueueReadBuffer(device->get_queue(), cpu.sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_back2, 0, 0, 0);
    BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
    logger.info() << cpu_back2;
    logger.info() << "Run kernel";
    int tmp = EVEN;
    //this is then force(Y_even, X_odd) == force(in1, in4)
    device->fermion_force_eo_device(cpu.in1_eo, cpu.in4_eo, device->get_gaugefield(), cpu.out, tmp, params.get_kappa() );
    
    logger.info() << "|force (even)|^2:";

    device->set_float_to_gaugemomentum_squarenorm_device(cpu.out, cpu.sqnorm);
    err = clEnqueueReadBuffer(device->get_queue(), cpu.sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_res, 0, 0, 0);
    BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
    logger.info() << cpu_res;
  } else {
    logger.info() << "|phi_odd_1|^2:";
    device->set_float_to_global_squarenorm_eoprec_device(cpu.in3, cpu.sqnorm);
    err = clEnqueueReadBuffer(device->get_queue(), cpu.sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_back3, 0, 0, 0);
    BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
    logger.info() << cpu_back3;
    logger.info() << "|phi_odd_2|^2:";
    device->set_float_to_global_squarenorm_eoprec_device(cpu.in4, cpu.sqnorm);
    err = clEnqueueReadBuffer(device->get_queue(), cpu.sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_back4, 0, 0, 0);
    BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
    logger.info() << cpu_back4;
    logger.info() << "Run kernel";
    
    int tmp = ODD;
    //this is then force(Y_odd, X_even) == force(in2, in3)
    device->fermion_force_eo_device(cpu.in2, cpu.in3, device->get_gaugefield(), cpu.out, tmp, params.get_kappa() );
   
    logger.info() << "|force (odd)|^2:";
    device->set_float_to_gaugemomentum_squarenorm_device(cpu.out, cpu.sqnorm);
    err = clEnqueueReadBuffer(device->get_queue(), cpu.sqnorm, CL_TRUE, 0, sizeof(hmc_float), &cpu_res, 0, 0, 0);
    BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
    logger.info() << cpu_res;
  }
  
  logger.info() << "Choosing reference value and acceptance precision";
  hmc_float ref_val = params.get_test_ref_value();
  logger.info() << "reference value:\t" << ref_val;
  hmc_float prec = params.get_solver_prec();
  logger.info() << "acceptance precision: " << prec;
  
  logger.info() << "Compare result to reference value";
  BOOST_REQUIRE_CLOSE(cpu_res, ref_val, prec);
  logger.info() << "Done";
  BOOST_MESSAGE("Test done");

}

BOOST_AUTO_TEST_SUITE( GF_UPDATE )

BOOST_AUTO_TEST_CASE( GF_UPDATE_1 ){
  test_gf_update("/gf_update_input_1");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_UPDATE )

BOOST_AUTO_TEST_CASE( F_UPDATE_1 ){
  test_f_update("/f_update_input_1");
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( F_GAUGE )

BOOST_AUTO_TEST_CASE( F_GAUGE_1 ){
  test_f_gauge("/f_gauge_input_1");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_GAUGE_TLSYM )

BOOST_AUTO_TEST_CASE( F_GAUGE_TLSYM_1 ){
  test_f_gauge_tlsym("/f_gauge_tlsym_input_1");
}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE( F_FERMION ) 

BOOST_AUTO_TEST_CASE( F_FERMION_1 ){
  test_f_fermion("/f_fermion_input_1");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( F_FERMION_EO ) 

BOOST_AUTO_TEST_CASE( F_FERMION_EO_1 ){
  //  test_f_fermion_eo("/f_fermion_eo_input_1");
}

BOOST_AUTO_TEST_SUITE_END()
