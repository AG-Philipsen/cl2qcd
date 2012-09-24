#include "../opencl_module_hmc.h"
#include "../gaugefield_hybrid.h"
#include "../meta/util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE gf_update
#include <boost/test/unit_test.hpp>

extern std::string const version;
std::string const version = "0.1";
std::string const exec_name = "gf_update_test";

class Device : public Opencl_Module_Hmc {

	cl_kernel testKernel;
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
	void get_gaugeobservables_from_task(int ntask, hmc_float * plaq, hmc_float * tplaq, hmc_float * splaq, hmc_complex * pol);
	void runTestKernel();

private:
	void fill_buffers();
	void clear_buffers();
	cl_mem in, out;
	cl_mem sqnorm;
	hmc_float * gm_in;
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

void fill_with_one(hmc_float * sf_in, int size)
{
	for(int i = 0; i < size; ++i) {
		sf_in[i] = 1.;
	}
	return;
}

void fill_with_random(hmc_float * sf_in, int size)
{
  prng_init(123456);
  for(int i = 0; i < size; ++i) {
    sf_in[i] = prng_double();
  }
  return;
}

void Dummyfield::fill_buffers()
{
	// don't invoke parent function as we don't require the original buffers
	cl_int err;
	cl_context context = opencl_modules[0]->get_context();
	int NUM_ELEMENTS_AE = meta::get_vol4d(get_parameters()) * NDIM * meta::get_su3algebrasize();
	gm_in = new hmc_float[NUM_ELEMENTS_AE];

	//use the variable use_cg to switch between cold and random input sf
	if(get_parameters().get_solver() == meta::Inputparameters::cg) {
		fill_with_one(gm_in, NUM_ELEMENTS_AE);
	} else {
		fill_with_random(gm_in, NUM_ELEMENTS_AE);
	}
	BOOST_REQUIRE(gm_in);

	Device * device = static_cast<Device*>(opencl_modules[0]);
	//create buffer for sf on device (and copy sf_in to both for convenience)

	in = clCreateBuffer(context, CL_MEM_READ_ONLY, device->get_gaugemomentum_buffer_size(), 0, &err);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	device->importGaugemomentumBuffer(in, reinterpret_cast<ae*>(gm_in));

	//create buffer for squarenorm on device
	sqnorm = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float), 0, &err);
}

void Device::fill_kernels()
{
	//one only needs some kernels up to now. to save time during compiling they are put in here by hand
	Opencl_Module_Hmc::fill_kernels();
}

void Dummyfield::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers
	clReleaseMemObject(in);
	clReleaseMemObject(sqnorm);

	delete[] gm_in;
}

void Device::clear_kernels()
{
	clReleaseKernel(testKernel);
	Opencl_Module::clear_kernels();
}

hmc_float Dummyfield::get_squarenorm()
{
	static_cast<Device*>(opencl_modules[0])->set_float_to_gaugemomentum_squarenorm_device(in, sqnorm);
	// get stuff from device
	hmc_float result;
	cl_int err = clEnqueueReadBuffer(*queue, sqnorm, CL_TRUE, 0, sizeof(hmc_float), &result, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	return result;
}

void Dummyfield::runTestKernel()
{
  hmc_float eps = 0.12;
  Device * device = static_cast<Device*>(opencl_modules[0]);
  static_cast<Device*>(opencl_modules[0])->md_update_gaugefield_device( in, device->get_gaugefield(), eps);
}

void Dummyfield::get_gaugeobservables_from_task(int ntask, hmc_float * plaq, hmc_float * tplaq, hmc_float * splaq, hmc_complex * pol)
{
	if( ntask < 0 || ntask > get_num_tasks() ) throw Print_Error_Message("devicetypes index out of range", __FILE__, __LINE__);
	opencl_modules[ntask]->gaugeobservables(plaq, tplaq, splaq, pol);
}

BOOST_AUTO_TEST_CASE( GF_UPDATE )
{
  logger.info() << "Test kernel";
  logger.info() << "\tmd_update_gaugefield";
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
  logger.info() << "|in|^2:";
  hmc_float cpu_back = cpu.get_squarenorm();
  logger.info() << cpu_back;
  cpu.runTestKernel();
  logger.info() << "gaugeobservables: ";
  cpu.print_gaugeobservables_from_task(0, 0);
  
  hmc_float plaq_cpu, tplaq_cpu, splaq_cpu;
  hmc_complex pol_cpu;
  cpu.get_gaugeobservables_from_task(0, &plaq_cpu, &tplaq_cpu, &splaq_cpu, &pol_cpu);
  BOOST_MESSAGE("Tested CPU");
  
  logger.info() << "Choosing reference value and acceptance precision";
  hmc_float ref_val = params.get_test_ref_value();
  logger.info() << "reference value:\t" << ref_val;
  hmc_float prec = params.get_solver_prec();  
  logger.info() << "acceptance precision: " << prec;

  logger.info() << "Compare result to reference value";
  BOOST_REQUIRE_CLOSE(plaq_cpu, ref_val, prec);
  logger.info() << "Done";
  BOOST_MESSAGE("Test done");
}
