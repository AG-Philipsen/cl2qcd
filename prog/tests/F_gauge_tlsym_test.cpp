#include "../opencl_module_hmc.h"
#include "../gaugefield_hybrid.h"
#include "../meta/util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Rectangles
#include <boost/test/unit_test.hpp>

extern std::string const version;
std::string const version = "0.1";
std::string const exec_name = "f_tlsym";

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

  hmc_float get_squarenorm();
  void runTestKernel();
  
private:
  void fill_buffers();
  void clear_buffers();
  cl_mem out;
  cl_mem sqnorm;
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

void fill_with_zero(hmc_float * sf_in, int size)
{
	for(int i = 0; i < size; ++i) {
		sf_in[i] = 0.;
	}
	return;
}

void Dummyfield::fill_buffers()
{
	// don't invoke parent function as we don't require the original buffers
	cl_int err;
	cl_context context = opencl_modules[0]->get_context();
	int NUM_ELEMENTS_AE = meta::get_vol4d(get_parameters()) * NDIM * meta::get_su3algebrasize();
	sf_out = new hmc_float[NUM_ELEMENTS_AE];
	fill_with_zero(sf_out, NUM_ELEMENTS_AE);

	Device * device = static_cast<Device*>(opencl_modules[0]);
	out = clCreateBuffer(context, CL_MEM_WRITE_ONLY, device->get_gaugemomentum_buffer_size(), 0, &err);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	device->importGaugemomentumBuffer(out, reinterpret_cast<ae*>(sf_out));

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
	clReleaseMemObject(out);
	clReleaseMemObject(sqnorm);
	delete[] sf_out;
}

void Device::clear_kernels()
{
	Opencl_Module::clear_kernels();
}

hmc_float Dummyfield::get_squarenorm()
{
	static_cast<Device*>(opencl_modules[0])->set_float_to_gaugemomentum_squarenorm_device(out, sqnorm);
	// get stuff from device
	hmc_float result;
	cl_int err = clEnqueueReadBuffer(*queue, sqnorm, CL_TRUE, 0, sizeof(hmc_float), &result, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	return result;
}

void Dummyfield::runTestKernel()
{
  Device * device = static_cast<Device*>(opencl_modules[0]);
  static_cast<Device*>(opencl_modules[0])->gauge_force_tlsym_device( device->get_gaugefield()  ,out);
}

BOOST_AUTO_TEST_CASE( F_GAUGE_TLSYM)
{
  logger.info() << "Test kernel";
  logger.info() << "\tf_tlsym";
  logger.info() << "against reference values";

  int param_expect = 4;
  logger.info() << "expect parameters:";
  logger.info() << "\texec_name\tinputfile\tgpu_usage\trec12_usage";
  //get number of parameters
  int num_par = boost::unit_test::framework::master_test_suite().argc;
  if(num_par < param_expect){
    logger.fatal() << "need more inputparameters! Got only " << num_par << ", expected " << param_expect << "! Aborting...";
    exit(-1);
  }

  //get input file that has been passed as an argument 
  const char*  inputfile =  boost::unit_test::framework::master_test_suite().argv[1];
  logger.info() << "inputfile used: " << inputfile;
  //get use_gpu = true/false that has been passed as an argument 
  const char*  gpu_opt =  boost::unit_test::framework::master_test_suite().argv[2];
  logger.info() << "GPU usage: " << gpu_opt;
  //get use_rec12 = true/false that has been passed as an argument 
  const char* rec12_opt =  boost::unit_test::framework::master_test_suite().argv[3];
  logger.info() << "rec12 usage: " << rec12_opt;

  logger.info() << "Init device";
  const char* _params_cpu[] = {"foo", inputfile, gpu_opt, rec12_opt};
  meta::Inputparameters params(param_expect, _params_cpu);
  Dummyfield cpu(params);
  logger.info() << "gaugeobservables: ";
  cpu.print_gaugeobservables_from_task(0, 0);
  cpu.runTestKernel();
  logger.info() << "|f_gauge|^2:";
  hmc_float cpu_res;
  cpu_res = cpu.get_squarenorm();
  logger.info() << cpu_res;

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

