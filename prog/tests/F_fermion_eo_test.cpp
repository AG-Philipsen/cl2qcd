#include "../opencl_module_hmc.h"
#include "../gaugefield_hybrid.h"

#include "../meta/util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE F_fermion_eo
#include <boost/test/unit_test.hpp>

extern std::string const version;
std::string const version = "0.1";
std::string const exec_name = "f_fermion_eo_test";

class Device : public Opencl_Module_Hmc {

	meta::Counter counter1, counter2, counter3, counter4;
public:
	Device(cl_command_queue queue, const meta::Inputparameters& params, int maxcomp, std::string double_ext, unsigned int dev_rank) : Opencl_Module_Hmc(params, &counter1, &counter2, &counter3, &counter4) {
		Opencl_Module_Hmc::init(queue, maxcomp, double_ext, dev_rank); /* init in body for proper this-pointer */
	};
	~Device() {
		finalize();
	};

	void fill_kernels();
	void clear_kernels();
};

class Dummyfield : public Gaugefield_hybrid {

public:	
  Dummyfield(meta::Inputparameters inputfile) : Gaugefield_hybrid(inputfile) {
    cl_device_type primary_device;
    switch ( inputfile.get_use_gpu() ) {
    case true :
      primary_device = CL_DEVICE_TYPE_GPU;
      break;
    case false :
      primary_device = CL_DEVICE_TYPE_CPU;
      break;
    }
    init(1, primary_device);
    meta::print_info_hmc(exec_name.c_str(), inputfile);
  };
  
	virtual void init_tasks();
	virtual void finalize_opencl();

	hmc_float get_squarenorm(int which);
	void runTestKernel(int evenodd);
	void reset_outfield();
private:
	void fill_buffers();
	void clear_buffers();
	cl_mem in1, in2, in3, in4, out;
	cl_mem sqnorm;
	spinor * sf_in1;
	spinor * sf_in2;
	spinor * sf_in3;
	spinor * sf_in4;
	hmc_float * sf_out;
};


void Dummyfield::init_tasks()
{
	opencl_modules = new Opencl_Module* [get_num_tasks()];
	opencl_modules[0] = new Device(queue[0], get_parameters(), get_max_compute_units(0), get_double_ext(0), 0);

	fill_buffers();
}

void Dummyfield::finalize_opencl()
{
	clear_buffers();
	Gaugefield_hybrid::finalize_opencl();
}

void fill_sf_with_one(spinor * sf_in1, int size)
{
	for(int i = 0; i < size; ++i) {
		sf_in1[i].e0.e0 = hmc_complex_one;
		sf_in1[i].e0.e1 = hmc_complex_one;
		sf_in1[i].e0.e2 = hmc_complex_one;
		sf_in1[i].e1.e0 = hmc_complex_one;
		sf_in1[i].e1.e1 = hmc_complex_one;
		sf_in1[i].e1.e2 = hmc_complex_one;
		sf_in1[i].e2.e0 = hmc_complex_one;
		sf_in1[i].e2.e1 = hmc_complex_one;
		sf_in1[i].e2.e2 = hmc_complex_one;
		sf_in1[i].e3.e0 = hmc_complex_one;
		sf_in1[i].e3.e1 = hmc_complex_one;
		sf_in1[i].e3.e2 = hmc_complex_one;
	}
	return;
}

void fill_sf_with_float(spinor * sf_in, int size, hmc_float val)
{
	for(int i = 0; i < size; ++i) {
		sf_in[i].e0.e0.re = val;
		sf_in[i].e0.e1.re = val;
		sf_in[i].e0.e2.re = val;
		sf_in[i].e1.e0.re = val;
		sf_in[i].e1.e1.re = val;
		sf_in[i].e1.e2.re = val;
		sf_in[i].e2.e0.re = val;
		sf_in[i].e2.e1.re = val;
		sf_in[i].e2.e2.re = val;
		sf_in[i].e3.e0.re = val;
		sf_in[i].e3.e1.re = val;
		sf_in[i].e3.e2.re = val;

		sf_in[i].e0.e0.im = val;
		sf_in[i].e0.e1.im = val;
		sf_in[i].e0.e2.im = val;
		sf_in[i].e1.e0.im = val;
		sf_in[i].e1.e1.im = val;
		sf_in[i].e1.e2.im = val;
		sf_in[i].e2.e0.im = val;
		sf_in[i].e2.e1.im = val;
		sf_in[i].e2.e2.im = val;
		sf_in[i].e3.e0.im = val;
		sf_in[i].e3.e1.im = val;
		sf_in[i].e3.e2.im = val;
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

void fill_with_zero(hmc_float * sf_in, int size)
{
	for(int i = 0; i < size; ++i) {
		sf_in[i] = 0.;
	}
	return;
}


void fill_sf_with_random(spinor * sf_in1, spinor * sf_in2, int size, int seed)
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

void Dummyfield::reset_outfield()
{
	static_cast<Device*>(opencl_modules[0])->importGaugemomentumBuffer(out, reinterpret_cast<ae*>(sf_out));
}

void Dummyfield::fill_buffers()
{
	// don't invoke parent function as we don't require the original buffers
	cl_int err;
	cl_context context = opencl_modules[0]->get_context();
	int NUM_ELEMENTS_SF;
	if(get_parameters().get_use_eo() == true) NUM_ELEMENTS_SF =  meta::get_eoprec_spinorfieldsize(get_parameters());
	else NUM_ELEMENTS_SF =  get_spinorfieldsize(get_parameters());

	int NUM_ELEMENTS_AE = meta::get_vol4d(get_parameters()) * NDIM * meta::get_su3algebrasize();

	sf_in1 = new spinor[NUM_ELEMENTS_SF];
	sf_in2 = new spinor[NUM_ELEMENTS_SF];
	sf_in3 = new spinor[NUM_ELEMENTS_SF];
	sf_in4 = new spinor[NUM_ELEMENTS_SF];
	sf_out = new hmc_float[NUM_ELEMENTS_AE];

	//use the variable use_cg to switch between cold and random input sf
	if(get_parameters().get_solver() == meta::Inputparameters::cg) {
		fill_sf_with_one(sf_in1, NUM_ELEMENTS_SF);
		fill_sf_with_one(sf_in2, NUM_ELEMENTS_SF);
		fill_sf_with_one(sf_in3, NUM_ELEMENTS_SF);
		fill_sf_with_one(sf_in4, NUM_ELEMENTS_SF);
	} else {
		fill_sf_with_random(sf_in1, sf_in2, NUM_ELEMENTS_SF, 123456);
		fill_sf_with_random(sf_in3, sf_in4, NUM_ELEMENTS_SF, 789101);
	}
	BOOST_REQUIRE(sf_in1);
	BOOST_REQUIRE(sf_in2);
	BOOST_REQUIRE(sf_in3);
	BOOST_REQUIRE(sf_in4);

	fill_with_zero(sf_out, NUM_ELEMENTS_AE);

	//create buffer for sf on device (and copy sf_in to both for convenience)

	Device * spinor_module = static_cast<Device*>(opencl_modules[0]);
	size_t sf_eoprec_buffer_size = spinor_module->get_eoprec_spinorfield_buffer_size();
	//create buffer for sf on device (and copy sf_in to both for convenience)

	in1 = clCreateBuffer(context, CL_MEM_READ_ONLY , sf_eoprec_buffer_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	in2 = clCreateBuffer(context, CL_MEM_READ_ONLY , sf_eoprec_buffer_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	in3 = clCreateBuffer(context, CL_MEM_READ_ONLY , sf_eoprec_buffer_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	in4 = clCreateBuffer(context, CL_MEM_READ_ONLY , sf_eoprec_buffer_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	spinor_module->copy_to_eoprec_spinorfield_buffer(in1, sf_in1);
	spinor_module->copy_to_eoprec_spinorfield_buffer(in2, sf_in2);
	spinor_module->copy_to_eoprec_spinorfield_buffer(in3, sf_in3);
	spinor_module->copy_to_eoprec_spinorfield_buffer(in4, sf_in4);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	out = clCreateBuffer(context, CL_MEM_WRITE_ONLY, spinor_module->get_gaugemomentum_buffer_size(), 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	this->reset_outfield();

	//create buffer for squarenorm on device
	sqnorm = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float), 0, &err);
}

void Device::fill_kernels()
{
	Opencl_Module_Hmc::fill_kernels();
}

void Dummyfield::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers
	clReleaseMemObject(in1);
	clReleaseMemObject(in2);
	clReleaseMemObject(in3);
	clReleaseMemObject(in4);
	clReleaseMemObject(out);
	clReleaseMemObject(sqnorm);

	delete[] sf_in1;
	delete[] sf_in2;
	delete[] sf_in3;
	delete[] sf_in4;
	delete[] sf_out;
}

void Device::clear_kernels()
{
	Opencl_Module_Hmc::clear_kernels();
}

hmc_float Dummyfield::get_squarenorm(int which)
{
	//which controlls if the in or out-vector is looked at
	if(which == 0) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_eoprec_device(in1, sqnorm);
	if(which == 1) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_eoprec_device(in2, sqnorm);
	if(which == 2) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_eoprec_device(in3, sqnorm);
	if(which == 3) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_eoprec_device(in4, sqnorm);
	if(which == 4) static_cast<Device*>(opencl_modules[0])->set_float_to_gaugemomentum_squarenorm_device(out, sqnorm);
	// get stuff from device
	hmc_float result;
	cl_int err = clEnqueueReadBuffer(*queue, sqnorm, CL_TRUE, 0, sizeof(hmc_float), &result, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	return result;
}

void Dummyfield::runTestKernel(int evenodd)
{
  //interprete Y = (in1, in2) X = (in3, in4)
  //Y_odd = in2, Y_even = in1, X_odd = in4, X_even = in3
  Device * device = static_cast<Device*>(opencl_modules[0]);
  if(evenodd == ODD) {
    //this is then force(Y_odd, X_even) == force(in2, in3)
    device->fermion_force_eo_device(in2, in3, device->get_gaugefield(), out, evenodd, get_parameters().get_kappa() );
  } else {
    //this is then force(Y_even, X_odd) == force(in1, in4)
    device->fermion_force_eo_device(in1, in4, device->get_gaugefield(), out, evenodd, get_parameters().get_kappa() );
  }
}

BOOST_AUTO_TEST_CASE( F_FERMION_EO )
{
  logger.info() << "Test kernel";
  logger.info() << "\tf_fermion_eo";
  logger.info() << "against reference value";

  //get input file that has been passed as an argument 
  const char* inputfile =  boost::unit_test::framework::master_test_suite().argv[1];
  logger.info() << "inputfile used: " << inputfile;
  //get use_gpu = true/false that has been passed as an argument 
  const char* gpu_opt =  boost::unit_test::framework::master_test_suite().argv[2];
  logger.info() << "GPU usage: " << gpu_opt;

  logger.info() << "Init device";
  const char* _params_cpu[] = {"foo", inputfile, gpu_opt};
  meta::Inputparameters params(3, _params_cpu);
  Dummyfield cpu(params);
  logger.info() << "gaugeobservables: ";
  cpu.print_gaugeobservables_from_task(0, 0);

  //switch according to "use_pointsource"
  hmc_float cpu_res;
  if(params.get_use_pointsource()){
    logger.info() << "|phi_even_1|^2:";
    hmc_float cpu_back = cpu.get_squarenorm(0);
    logger.info() << cpu_back;
    logger.info() << "|phi_even_2|^2:";
    hmc_float cpu_back2 = cpu.get_squarenorm(1);
    logger.info() << cpu_back2;
    cpu.runTestKernel(EVEN);
    logger.info() << "|force (even)|^2:";
    cpu_res = cpu.get_squarenorm(4);    
    logger.info() << cpu_res;
  } else {
    logger.info() << "|phi_odd_1|^2:";
    hmc_float cpu_back3 = cpu.get_squarenorm(2);
    logger.info() << cpu_back3;
    logger.info() << "|phi_odd_2|^2:";
    hmc_float cpu_back4 = cpu.get_squarenorm(3);
    logger.info() << cpu_back4;
    logger.info() << "Run kernel";
    cpu.runTestKernel(ODD);
    logger.info() << "|force (odd)|^2:";
    cpu_res = cpu.get_squarenorm(4);
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
