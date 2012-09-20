#include "../opencl_module_spinors.h"
#include "../gaugefield_hybrid.h"

#include "../meta/util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Memory Access using Complex
#include <boost/test/unit_test.hpp>

extern std::string const version;
std::string const version = "0.1";

class Device : public Opencl_Module {

	cl_kernel fillComplex;
	cl_kernel readComplex;

public:
	Device(cl_command_queue queue, const meta::Inputparameters& params, hardware::Device * device, int maxcomp, std::string double_ext, unsigned int dev_rank)
		: Opencl_Module(params, device) {
		Opencl_Module::init(queue, maxcomp, double_ext, dev_rank); /* init in body for proper this-pointer */
	};
	~Device() {
		finalize();
	};

	void runFillKernel(cl_mem out, hmc_complex value);
	void runReadKernel(cl_mem out, cl_mem in);
	void fill_kernels();
	void clear_kernels();
};

class Dummyfield : public Gaugefield_hybrid {

public:
	Dummyfield(cl_device_type device_type, const hardware::System * system)
		: Gaugefield_hybrid(system) {
		init(1, device_type);
	};

	virtual void init_tasks();
	virtual void finalize_opencl();

	void verifyFill(const hmc_complex value);
	void runFillKernel(const hmc_complex value);
	void verifyRead();
	void runReadKernel();

private:
	void verify(cl_float2, cl_float2);
	void verify(hmc_complex, hmc_complex);
	void fill_buffers();
	void clear_buffers();
	hmc_complex * h_complex;
	cl_mem d_complex;
	cl_mem d_readComplex;
	cl_mem d_float2;
};

BOOST_AUTO_TEST_CASE( CPU )
{
	const char* _params_cpu[] = {"foo", "--use_gpu=false"};
	meta::Inputparameters params_cpu(2, _params_cpu);
	hardware::System system_cpu(params_cpu);
	Dummyfield dummy(CL_DEVICE_TYPE_CPU, &system_cpu);
	dummy.runFillKernel(hmc_complex_zero);
	dummy.verifyFill(hmc_complex_zero);
	dummy.runFillKernel(hmc_complex_one);
	dummy.verifyFill(hmc_complex_one);
	BOOST_MESSAGE("Tested CPU");
}

BOOST_AUTO_TEST_CASE( GPU )
{
	const char* _params_gpu[] = {"foo", "--use_gpu=true"};
	meta::Inputparameters params_gpu(2, _params_gpu);
	hardware::System system_gpu(params_gpu);
	Dummyfield dummy(CL_DEVICE_TYPE_GPU, &system_gpu);
	dummy.runFillKernel(hmc_complex_zero);
	dummy.verifyFill(hmc_complex_zero);
	dummy.runFillKernel(hmc_complex_one);
	dummy.verifyFill(hmc_complex_one);
	BOOST_MESSAGE("Tested GPU");
}

void Dummyfield::init_tasks()
{
	opencl_modules = new Opencl_Module* [get_num_tasks()];
	opencl_modules[0] = new Device(queue[0], get_parameters(), get_device_for_task(0), get_max_compute_units(0), get_double_ext(0), 0);

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

	h_complex = new hmc_complex[meta::get_vol4d(get_parameters())];
	BOOST_REQUIRE(h_complex);

	for(int i = 0; i < meta::get_vol4d(get_parameters()); ++i) {
		hmc_complex tmp = { (hmc_float) i, (hmc_float) (meta::get_vol4d(get_parameters()) - i) };
		h_complex[i] = tmp;
	}
	d_readComplex = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, meta::get_vol4d(get_parameters()) * sizeof(hmc_complex), h_complex, &err);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	d_float2 = clCreateBuffer(context, CL_MEM_WRITE_ONLY, meta::get_vol4d(get_parameters()) * sizeof(cl_float2), 0, &err);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	d_complex = clCreateBuffer(context, CL_MEM_WRITE_ONLY, meta::get_vol4d(get_parameters()) * sizeof(hmc_complex), 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
}

void Device::fill_kernels()
{
	Opencl_Module::fill_kernels();

	fillComplex = createKernel("fillComplex") << basic_opencl_code << "tests/complex_mem_access.cl";
	readComplex = createKernel("readComplex") << basic_opencl_code << "tests/complex_mem_access.cl";
}

void Dummyfield::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	clReleaseMemObject(d_readComplex);
	clReleaseMemObject(d_float2);
	clReleaseMemObject(d_complex);

	delete[] h_complex;
}

void Device::clear_kernels()
{
	clReleaseKernel(fillComplex);
	Opencl_Module::clear_kernels();
}

void Device::runFillKernel(cl_mem out, hmc_complex value)
{
	cl_int err;
	err = clSetKernelArg(fillComplex, 0, sizeof(cl_mem), &out);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(fillComplex, 1, sizeof(hmc_complex), &value);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	enqueueKernel(fillComplex, 1024);
}

void Device::runReadKernel(cl_mem out, cl_mem in)
{
	cl_int err;
	err = clSetKernelArg(readComplex, 0, sizeof(cl_mem), &out);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(readComplex, 1, sizeof(cl_mem), &in);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	enqueueKernel(readComplex, 1024);
}

void Dummyfield::verify(hmc_complex left, hmc_complex right)
{
	BOOST_REQUIRE_EQUAL(left.re, right.re);
	BOOST_REQUIRE_EQUAL(left.im, right.im);
}


void Dummyfield::verifyFill(const hmc_complex value)
{
	// get stuff from device
	cl_int err = clEnqueueReadBuffer(*queue, d_complex, CL_TRUE, 0, meta::get_vol4d(get_parameters()) * sizeof(hmc_complex), h_complex, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);

	for(int i = 0; i < meta::get_vol4d(get_parameters()); ++i) {
		verify(value, h_complex[i]);
	}
}

void Dummyfield::runFillKernel(const hmc_complex value)
{
	static_cast<Device*>(opencl_modules[0])->runFillKernel(d_complex, value);
}

// usually one shouldn't verify flaots in this way, but as we actually only store ints ...
void Dummyfield::verify(cl_float2 left, cl_float2 right)
{
	BOOST_REQUIRE_EQUAL(left.s[0], right.s[0]);
	BOOST_REQUIRE_EQUAL(left.s[1], right.s[1]);
}

void Dummyfield::verifyRead()
{
	// get stuff from device
	cl_float2 * h_float2 = new cl_float2[meta::get_vol4d(get_parameters())];
	cl_int err = clEnqueueReadBuffer(*queue, d_float2, CL_TRUE, 0, meta::get_vol4d(get_parameters()) * sizeof(cl_float2), h_float2, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);

	for(int i = 0; i < meta::get_vol4d(get_parameters()); ++i) {
		cl_float2 ref = {(cl_float) i, (cl_float) (meta::get_vol4d(get_parameters()) - i)};
		verify(ref, h_float2[i]);
	}

	delete[] h_float2;
}

void Dummyfield::runReadKernel()
{
	static_cast<Device*>(opencl_modules[0])->runReadKernel(d_float2, d_readComplex);
}

