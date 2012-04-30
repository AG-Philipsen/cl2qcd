#include "../opencl_module_fermions.h"
#include "../gaugefield_hybrid.h"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Fermionmatrix
#include <boost/test/unit_test.hpp>

extern std::string const version;
std::string const version = "0.1";

class Device : public Opencl_Module {

	cl_kernel testKernel;

public:
	Device(cl_command_queue queue, inputparameters* params, int maxcomp, string double_ext, unsigned int dev_rank) : Opencl_Module() {
		Opencl_Module::init(queue, params, maxcomp, double_ext, dev_rank); /* init in body for proper this-pointer */
	};
	~Device() {
		finalize();
	};

	void runTestKernel(int gs, int ls);
	void fill_kernels();
	void clear_kernels();
};

class Dummyfield : public Gaugefield_hybrid {

public:
	Dummyfield(cl_device_type device_type) : Gaugefield_hybrid() {

		init(1, device_type, &params);
	};

	virtual void init_tasks();
	virtual void finalize_opencl();

	void runTestKernel();

private:
	void fill_buffers();
	void clear_buffers();
	inputparameters params;
};


BOOST_AUTO_TEST_CASE( GEOMETRY )
{
	logger.info() << "Init CPU device";
	Dummyfield dummy(CL_DEVICE_TYPE_CPU);
	dummy.runTestKernel();
	logger.info() << "Init GPU device";
	Dummyfield dummy2(CL_DEVICE_TYPE_GPU);
	dummy2.runTestKernel();
	BOOST_MESSAGE("Tested GEOMETRY");
}

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
	return;
}

void Device::fill_kernels()
{
	//one only needs some kernels up to now. to save time during compiling they are put in here by hand
	Opencl_Module::fill_kernels();

	//to this end, one has to set the needed files by hand
	basic_opencl_code = ClSourcePackage() << "opencl_header.cl" << "operations_geometry.cl";
	testKernel = createKernel("geometry_test") << basic_opencl_code << "tests/geometry_test.cl";
}

void Dummyfield::clear_buffers()
{
}

void Device::clear_kernels()
{
	clReleaseKernel(testKernel);
}

void Device::runTestKernel(int gs, int ls)
{
	enqueueKernel(testKernel, gs, ls);
}

void Dummyfield::runTestKernel()
{
	int gs, ls;
	if(opencl_modules[0]->get_device_type() == CL_DEVICE_TYPE_GPU) {
		gs = get_parameters()->get_eoprec_spinorfieldsize();
		ls = 64;
	} else if(opencl_modules[0]->get_device_type() == CL_DEVICE_TYPE_CPU) {
		gs = opencl_modules[0]->get_max_compute_units();
		ls = 1;
	}
	logger.info() << "test kernel with global_work_size: " << gs << " and local_work_size: " << ls;
	static_cast<Device*>(opencl_modules[0])->runTestKernel(gs, ls);
}

