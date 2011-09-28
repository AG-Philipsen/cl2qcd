#include "../opencl_module.h"
#include "../gaugefield_hybrid.h"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SU2 SU3 Extend
#include <boost/test/unit_test.hpp>

Random rnd(15);
extern std::string const version;
std::string const version = "0.1";

class Device : public Opencl_Module {

	cl_kernel testKernel_1;
	cl_kernel testKernel_2;

public:
	Device(cl_command_queue queue, inputparameters* params, int maxcomp, string double_ext) : Opencl_Module() {
		Opencl_Module::init(queue, 0, params, maxcomp, double_ext); /* init in body for proper this-pointer */
	};
	~Device() {
		finalize();
	};

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

private:
	inputparameters params;
};

BOOST_AUTO_TEST_CASE( GPU )
{
	Dummyfield dummy(CL_DEVICE_TYPE_GPU);
	BOOST_MESSAGE("Tested GPU");
}

void Dummyfield::init_tasks()
{
	opencl_modules = new Opencl_Module* [get_num_tasks()];
	opencl_modules[0] = new Device(queue[0], get_parameters(), get_max_compute_units(0), get_double_ext(0));
}

void Dummyfield::finalize_opencl()
{
	Gaugefield_hybrid::finalize_opencl();
}

void Device::fill_kernels()
{
	testKernel_1 = createKernel("saxsbypz_1") << "tests/saxsbypz.cl";
	testKernel_2 = createKernel("saxsbypz_2") << "tests/saxsbypz.cl";
}

void Device::clear_kernels()
{
	clReleaseKernel(testKernel_1);
	clReleaseKernel(testKernel_2);
}

