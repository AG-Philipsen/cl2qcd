#include "../opencl_module_hmc.h"
#include "../gaugefield_hybrid.h"
#include "../meta/util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Rectangles
#include <boost/test/unit_test.hpp>

extern std::string const version;
std::string const version = "0.1";

class TestDevice : public Opencl_Module {
public:
  //  TestDevice(cl_command_queue queue, const meta::Inputparameters& params, int maxcomp, std::string double_ext, unsigned int dev_rank) : Opencl_Module(params) {
  TestDevice(const meta::Inputparameters& params, hardware::Device * device) : Opencl_Module(params, device) {
    Opencl_Module::init(); /* init in body for proper this-pointer */
  };
  ~TestDevice() {
    finalize();
  };
  void fill_kernels();
  void clear_kernels();
};

void TestDevice::fill_kernels()
{
	Opencl_Module::fill_kernels();
}

void TestDevice::clear_kernels()
{
	Opencl_Module::clear_kernels();
}


class TestGaugefield : public Gaugefield_hybrid {
public:
	TestGaugefield(const hardware::System * system) : Gaugefield_hybrid(system) {
		auto inputfile = system->get_inputparameters();
		init(1, inputfile.get_use_gpu() ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU);
		std::string name = "test program";
		meta::print_info_hmc(name.c_str(), inputfile);
		logger.info() << "gaugeobservables: ";
		this->print_gaugeobservables_from_task(0, 0);
	};
	virtual void init_tasks() override;
	virtual void finalize_opencl() override;

	TestDevice * get_device();

private:
	void fill_buffers();
	void clear_buffers();
	cl_mem rect_value;
};

void TestGaugefield::init_tasks()
{
	opencl_modules = new Opencl_Module* [get_num_tasks()];
	opencl_modules[0] = new TestDevice(get_parameters(), get_device_for_task(0));

	fill_buffers();
}

TestDevice* TestGaugefield::get_device()
{
	return static_cast<TestDevice*>(opencl_modules[0]);
}

void TestGaugefield::finalize_opencl()
{
	clear_buffers();
	Gaugefield_hybrid::finalize_opencl();
}

void TestGaugefield::fill_buffers()
{
	// don't invoke parent function as we don't require the original buffers
	cl_int err;
	cl_context context = opencl_modules[0]->get_context();
	rect_value = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float), 0, &err);
}

void TestGaugefield::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers
	clReleaseMemObject(rect_value);
}

meta::Inputparameters create_parameters()
{
  const int param_expect = 4;
  logger.info() << "expect parameters:";
  logger.info() << "\texec_name\tinputfile\tgpu_usage\trec12_usage";
  std::string inputfile, gpu_opt , rec12_opt;
  //get number of parameters                                                                                                                                                                                                                 
  int num_par = boost::unit_test::framework::master_test_suite().argc;
  if(num_par < param_expect){
    logger.fatal() << "Got only " << num_par << " Inputparameters, expected " << param_expect << "! Use inputfile values instead!";
  } else if(num_par > param_expect){
    logger.warn() << "Got " << num_par << " Inputparameters, expected only " << param_expect << "! Use only the first " << param_expect << " values";
  }
  
  switch(num_par){
  case 0:
    logger.fatal() << "Something went terribly wrong! Did not even get executable name! Aborting...";
    exit(-1);
    break;
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
    logger.info() << rec12_opt;
    break;
  default:
    //get input file that has been passed as an argument                                                                                                                                                                                     
    inputfile =  boost::unit_test::framework::master_test_suite().argv[1];
    logger.info() << "inputfile used: " << inputfile;
    //get use_gpu = true/false that has been passed as an argument                                                                                                                                                                           
    gpu_opt =  boost::unit_test::framework::master_test_suite().argv[2];
    logger.info() << "GPU usage: " << gpu_opt;
    //get use_rec12 = true/false that has been passed as an argument                                                                                                                                                                         
    rec12_opt =  boost::unit_test::framework::master_test_suite().argv[3];
    logger.info() << "rec12 usage: " << rec12_opt;
    break;
  }
  const char* _params_cpu[] = {"foo", inputfile.c_str(), gpu_opt.c_str() , rec12_opt.c_str()};
  meta::Inputparameters params(num_par, _params_cpu);
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
	TestGaugefield cpu(&system);

	logger.info() << "calc rectangles value:";
	hmc_float cpu_rect;
	TestDevice * device = cpu.get_device();
	device->gaugeobservables_rectangles(device->get_gaugefield(), &cpu_rect);
	logger.info() << cpu_rect;
	logger.info() << "Finalize device";
	cpu.finalize();

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
