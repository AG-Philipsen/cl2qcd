#include "../gaugefield_hybrid.h"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SU2 SU3 Extend
#include <boost/test/unit_test.hpp>

extern std::string const version;
std::string const version = "0.1";

class Device : public Opencl_Module {

	cl_kernel testKernel;
	void fill_kernels();
	void clear_kernels();

protected:
	virtual size_t get_read_write_size(const std::string&) const {
		return 0;
	};
	virtual uint64_t get_flop_size(const std::string&) const {
		return 0;
	};
public:
	Device(const meta::Inputparameters& params, hardware::Device * device) : Opencl_Module(params, device) {
		fill_kernels();
	};
	~Device() {
		clear_kernels();
	};

};

const meta::Inputparameters INPUT(0, 0);

class Dummyfield : public Gaugefield_hybrid {

public:
	Dummyfield(cl_device_type device_type, const hardware::System * system)
		: Gaugefield_hybrid(system) {
		init(1, device_type);
	};

	virtual void init_tasks();
	virtual void finalize_opencl();
};

BOOST_AUTO_TEST_CASE( GPU )
{
	hardware::System system(INPUT);
	Dummyfield dummy(CL_DEVICE_TYPE_GPU, &system);
	BOOST_MESSAGE("Tested GPU");
}

void Dummyfield::init_tasks()
{
	opencl_modules = new Opencl_Module* [get_num_tasks()];
	opencl_modules[0] = new Device(get_parameters(), get_device_for_task(0));
}

void Dummyfield::finalize_opencl()
{
	Gaugefield_hybrid::finalize_opencl();
}

void Device::fill_kernels()
{
#ifdef _USEDOUBLEPREC_
	logger.info() << "init saxsbypz kernels to see if the compiler does strange things...";
	logger.info() << "\thmc_float = double";
	testKernel = createKernel("saxsbypz_1") << "tests/saxsbypz.cl";
	testKernel = createKernel("saxsbypz_2") << "tests/saxsbypz.cl";
#else
	logger.info() << "\thmc_float = single";
	testKernel = createKernel("saxsbypz_1") << "tests/saxsbypz_float.cl";
	testKernel = createKernel("saxsbypz_2") << "tests/saxsbypz_float.cl";
#endif
}

void Device::clear_kernels()
{
	clReleaseKernel(testKernel);
}

