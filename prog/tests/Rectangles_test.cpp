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
	Device(cl_command_queue queue, const meta::Inputparameters& params, int maxcomp, std::string double_ext, unsigned int dev_rank) : Opencl_Module_Hmc(params, &counter1, &counter2, &counter3, &counter4) {
		Opencl_Module_Hmc::init(queue, maxcomp, double_ext, dev_rank); /* init in body for proper this-pointer */
	};
	~Device() {
		finalize();
	};
	void fill_kernels();
	void clear_kernels();
};
/*
const std::string SOURCEFILE = std::string(SOURCEDIR)
#ifdef _USEDOUBLEPREC_
                               + "/tests/f_rectangles_input_1";
#else
                               + "/tests/f_rectangles_input_1_single";
#endif
const char * PARAMS[] = {"foo", SOURCEFILE.c_str()};
const meta::Inputparameters INPUT(2, PARAMS);

const std::string SOURCEFILE_REC12 = std::string(SOURCEDIR)
#ifdef _USEDOUBLEPREC_
                               + "/tests/f_rectangles_input_rec12";
#else
                               + "/tests/f_rectangles_input_rec12_single";
#endif
const char * PARAMS_REC12[] = {"foo", SOURCEFILE_REC12.c_str()};
const meta::Inputparameters INPUT_REC12(2, PARAMS_REC12);
*/


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

  hmc_float get_rect();
  
private:
  void fill_buffers();
  void clear_buffers();
  cl_mem rect_value;
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

BOOST_AUTO_TEST_CASE( RECTANGLES )
{
  logger.info() << "Test kernel";
  logger.info() << "\trectangles";
  logger.info() << "against reference value";

  //get input file that has been passed as an argument 
  const char* inputfile =  boost::unit_test::framework::master_test_suite().argv[1];
  logger.info() << "inputfile used: " << inputfile;
  //get use_gpu = true/false that has been passed as an argument 
  const char* gpu_opt =  boost::unit_test::framework::master_test_suite().argv[2];
  logger.info() << "GPU usage: " << gpu_opt;

  logger.info() << "Init device";
  const char* _params[] = {"foo", inputfile, gpu_opt};
  meta::Inputparameters params(3, _params);
  Dummyfield cpu(params);
  logger.info() << "gaugeobservables: ";
  cpu.print_gaugeobservables_from_task(0, 0);
  logger.info() << "calc rectangles value:";
  hmc_float cpu_rect = cpu.get_rect();
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

