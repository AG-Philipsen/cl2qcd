#include "../opencl_module_hmc.h"
#include "../gaugefield_hybrid.h"

#include "../meta/util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE staple_test
#include <boost/test/unit_test.hpp>

extern std::string const version;
std::string const version = "0.1";
std::string const exec_name = "staple_test";

class Device : public Opencl_Module {

	cl_kernel testKernel;
public:
	Device(cl_command_queue queue, const meta::Inputparameters& params, int maxcomp, std::string double_ext, unsigned int dev_rank) : Opencl_Module(params) {
		Opencl_Module::init(queue, maxcomp, double_ext, dev_rank); /* init in body for proper this-pointer */
	};
	~Device() {
		finalize();
	};

	void runTestKernel(cl_mem gf, cl_mem out, int gs, int ls);
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
  hmc_float runTestKernel();
  
private:
  void fill_buffers();
  void clear_buffers();
  cl_mem out;
  hmc_float * host_out;
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

void Dummyfield::fill_buffers()
{
	// don't invoke parent function as we don't require the original buffers
	cl_int err;
	cl_context context = opencl_modules[0]->get_context();
	int NUM_ELEMENTS = meta::get_vol4d(get_parameters());
	host_out = new hmc_float[NUM_ELEMENTS];
	BOOST_REQUIRE(host_out);

	size_t buf_size = NUM_ELEMENTS * sizeof(hmc_float);
	out = clCreateBuffer(context, CL_MEM_READ_ONLY , buf_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
}

void Device::fill_kernels()
{
	Opencl_Module::fill_kernels();
	testKernel = createKernel("staple_test") << basic_opencl_code  << "/tests/staple_test.cl";
}

void Dummyfield::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers
	clReleaseMemObject(out);
	delete[] host_out;
}

void Device::clear_kernels()
{
	clReleaseKernel(testKernel);
	Opencl_Module::clear_kernels();
}

void Device::runTestKernel(cl_mem gf, cl_mem out, int gs, int ls)
{
	cl_int err;
	err = clSetKernelArg(testKernel, 0, sizeof(cl_mem), &gf);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 1, sizeof(cl_mem), &out );
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);

	enqueueKernel(testKernel, gs, ls);
}

hmc_float Dummyfield::runTestKernel()
{
	hmc_float res = 0;
	int gs, ls;
	if(opencl_modules[0]->get_device_type() == CL_DEVICE_TYPE_GPU) {
		gs = meta::get_vol4d(get_parameters());
		ls = 64;
	} else {
		gs = opencl_modules[0]->get_max_compute_units();
		ls = 1;
	}
	Device * device = static_cast<Device*>(opencl_modules[0]);
	device->runTestKernel(device->get_gaugefield(), out, gs, ls);

	int NUM_ELEMENTS = meta::get_vol4d(get_parameters());
	//copy the result of the kernel to host
	size_t size = NUM_ELEMENTS * sizeof(hmc_float);
	cl_int clerr = clEnqueueReadBuffer(opencl_modules[0]->get_queue(), out, CL_TRUE, 0, size, host_out, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueReadBuffer", __FILE__, __LINE__);

	//sum up all elements in the result buffer
	for(int i = 0; i < NUM_ELEMENTS; i++) {
		res += host_out[i];
	}
	return res;
}


BOOST_AUTO_TEST_CASE( STAPLE_TEST )
{
  logger.info() << "Test kernel";
  logger.info() << "\tcalc_staple";
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
  logger.info() << "Run kernel";
  logger.info() << "running test kernel";
  hmc_float cpu_res = cpu.runTestKernel();
  logger.info() << "result:";
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

