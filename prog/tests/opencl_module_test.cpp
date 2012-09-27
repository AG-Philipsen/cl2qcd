#include "../opencl_module_hmc.h"
#include "../gaugefield_hybrid.h"
#include "../meta/util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE
#include <boost/test/unit_test.hpp>

//some functionality
#include "test_util.h"

extern std::string const version;
std::string const version = "0.1";

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

  Opencl_Module * get_device();

private:
	void fill_buffers();
	void clear_buffers();
	cl_mem rect_value;
};

void TestGaugefield::init_tasks()
{
	opencl_modules = new Opencl_Module* [get_num_tasks()];
	//here we want to test Opencl_Module
	opencl_modules[0] = new Opencl_Module(get_parameters(), get_device_for_task(0));
	opencl_modules[0]->init();

	fill_buffers();
}

Opencl_Module* TestGaugefield::get_device()
{
	return static_cast<Opencl_Module*>(opencl_modules[0]);
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

void test_rectangles(std::string inputfile)
{
  std::string kernelName = "rectangles";
  printKernelInfo(kernelName);
  logger.info() << "Init device";
  meta::Inputparameters params = create_parameters(inputfile);
  hardware::System system(params);
  TestGaugefield cpu(&system);

  logger.info() << "calc rectangles value:";
  hmc_float cpu_rect;
  Opencl_Module * device = cpu.get_device();
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
  BOOST_MESSAGE("Test done");
}

BOOST_AUTO_TEST_SUITE ( RECTANGLES )

BOOST_AUTO_TEST_CASE( RECTANGLES_1 )
{
  test_rectangles("/rectangles_input_1");
}

BOOST_AUTO_TEST_CASE( RECTANGLES_2 )
{
  test_rectangles("/rectangles_input_2");
}

BOOST_AUTO_TEST_SUITE_END()
