#include "../opencl_module_hmc.h"
#include "../gaugefield_hybrid.h"
#include "../meta/util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Rectangles
#include <boost/test/unit_test.hpp>

extern std::string const version;
std::string const version = "0.1";
std::string const exec_name = "rectangles";

class Device : public Opencl_Module_Hmc {
	meta::Counter counter1, counter2, counter3, counter4;
public:
	Device(const meta::Inputparameters& params, hardware::Device * device) : Opencl_Module_Hmc(params, device, &counter1, &counter2, &counter3, &counter4) {
		Opencl_Module_Hmc::init(); /* init in body for proper this-pointer */
	};
	~Device() {
		finalize();
	};
	void fill_kernels();
	void clear_kernels();
};

class Dummyfield : public Gaugefield_hybrid {
public:
	Dummyfield(const hardware::System * system) : Gaugefield_hybrid(system) {
		auto inputfile = system->get_inputparameters();
		init(1, inputfile.get_use_gpu() ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU);
		meta::print_info_hmc(exec_name.c_str(), inputfile);
		logger.info() << "gaugeobservables: ";
		this->print_gaugeobservables_from_task(0, 0);
	};
	virtual void init_tasks() override;
	virtual void finalize_opencl() override;

	hmc_float get_rect();
	Device * get_device();

private:
	void fill_buffers();
	void clear_buffers();
	cl_mem rect_value;
};

void Dummyfield::init_tasks()
{
	opencl_modules = new Opencl_Module* [get_num_tasks()];
	opencl_modules[0] = new Device(get_parameters(), get_device_for_task(0));

	fill_buffers();
}

Device* Dummyfield::get_device(){
  return static_cast<Device*>(opencl_modules[0]);
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
	rect_value = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float), 0, &err);
}

void Device::fill_kernels()
{
	Opencl_Module_Hmc::fill_kernels();
}

void Dummyfield::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers
	clReleaseMemObject(rect_value);
}

void Device::clear_kernels()
{
	Opencl_Module::clear_kernels();
}

hmc_float Dummyfield::get_rect()
{
	hmc_float rect_out;
	Device * device = static_cast<Device*>(opencl_modules[0]);
	device->gaugeobservables_rectangles(device->get_gaugefield(), &rect_out);
	return rect_out;
}

meta::Inputparameters create_parameters()
{
  int param_expect = 4;
  logger.info() << "expect parameters:";
  logger.info() << "\texec_name\tinputfile\tgpu_usage\trec12_usage";
  //get number of parameters
  const char* inputfile;
  const char* gpu_opt;
  const char* rec12_opt;

  int num_par = boost::unit_test::framework::master_test_suite().argc;
  if(num_par < param_expect){
    logger.fatal() << "need more inputparameters! Got only " << num_par << ", expected " << param_expect << "! Use inputfile values instead!";
    //exit(-1);
    switch(num_par){
    case 1:
      logger.fatal() << "Need at least an inputfile! Aborting...";
      exit(-1);
      break;
    case 2:
      //get input file that has been passed as an argument 
      inputfile =  boost::unit_test::framework::master_test_suite().argv[1];
      logger.info() << "inputfile used: " << inputfile;
      gpu_opt = "";
      rec12_opt = "";
      break;
    case 3:
      //get input file that has been passed as an argument 
      inputfile =  boost::unit_test::framework::master_test_suite().argv[1];
      logger.info() << "inputfile used: " << inputfile;
      //get use_gpu = true/false that has been passed as an argument 
      gpu_opt =  boost::unit_test::framework::master_test_suite().argv[2];
      logger.info() << "GPU usage: " << gpu_opt;
      rec12_opt = "";
      break;
    }
  }
  else{
    //get input file that has been passed as an argument 
    inputfile =  boost::unit_test::framework::master_test_suite().argv[1];
    logger.info() << "inputfile used: " << inputfile;
    //get use_gpu = true/false that has been passed as an argument 
    gpu_opt =  boost::unit_test::framework::master_test_suite().argv[2];
    logger.info() << "GPU usage: " << gpu_opt;
    //get use_rec12 = true/false that has been passed as an argument 
    rec12_opt =  boost::unit_test::framework::master_test_suite().argv[3];
    logger.info() << "rec12 usage: " << rec12_opt;
  }
  const char* _params_cpu[] = {"foo", inputfile, gpu_opt, rec12_opt};
  meta::Inputparameters params(param_expect, _params_cpu);
  return params;
}


BOOST_AUTO_TEST_CASE( RECTANGLES )
{
  logger.info() << "Test kernel";
  logger.info() << "\trectangles";
  logger.info() << "against reference value";

  logger.info() << "Init device";
  meta::Inputparameters params = create_parameters();
  hardware::System system(params);
  Dummyfield cpu(&system);

  logger.info() << "calc rectangles value:";
  hmc_float cpu_rect;
  Device * device = cpu.get_device();
  device->gaugeobservables_rectangles(device->get_gaugefield(), &cpu_rect);
  logger.info() << cpu_rect;

  logger.info() << "Choosing reference value and acceptance precision";
  hmc_float ref_val = params.get_test_ref_value();
  logger.info() << "reference value:\t" << ref_val;
  hmc_float prec = params.get_solver_prec();
  logger.info() << "acceptance precision: " << prec;

  logger.info() << "Compare result to reference value";
  BOOST_REQUIRE_CLOSE(cpu_rect, ref_val, prec);
  logger.info() << "Done";
  BOOST_MESSAGE("Test done");
}
