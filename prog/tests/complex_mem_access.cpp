#include "../opencl_module_spinors.h"
#include "../gaugefield_hybrid.h"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Memory Access using Complex
#include <boost/test/unit_test.hpp>

Random rnd(15);
extern std::string const version;
std::string const version = "0.1";

class Device : public Opencl_Module {

	cl_kernel fillComplex;
	cl_kernel readComplex;

public:
	Device(cl_command_queue queue, inputparameters* params, int maxcomp, string double_ext) : Opencl_Module() {
		Opencl_Module::init(queue, 0, params, maxcomp, double_ext); /* init in body for proper this-pointer */
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
	Dummyfield(cl_device_type device_type) : Gaugefield_hybrid() {
		init(1, device_type, &params);
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
	inputparameters params;
	hmc_complex * h_complex;
	cl_mem d_complex;
	cl_mem d_readComplex;
	cl_mem d_float2;
};

BOOST_AUTO_TEST_CASE( CPU )
{
	Dummyfield dummy(CL_DEVICE_TYPE_CPU);
	dummy.runFillKernel(hmc_complex_zero);
	dummy.verifyFill(hmc_complex_zero);
	dummy.runFillKernel(hmc_complex_one);
	dummy.verifyFill(hmc_complex_one);
	BOOST_MESSAGE("Tested CPU");
}

BOOST_AUTO_TEST_CASE( GPU )
{
	Dummyfield dummy(CL_DEVICE_TYPE_GPU);
	dummy.runFillKernel(hmc_complex_zero);
	dummy.verifyFill(hmc_complex_zero);
	dummy.runFillKernel(hmc_complex_one);
	dummy.verifyFill(hmc_complex_one);
	BOOST_MESSAGE("Tested GPU");
}

void Dummyfield::init_tasks()
{
	opencl_modules = new Opencl_Module* [get_num_tasks()];
	opencl_modules[0] = new Device(queue[0], get_parameters(), get_max_compute_units(0), get_double_ext(0));

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

	h_complex = new hmc_complex[params.get_vol4d()];
	BOOST_REQUIRE(h_complex);

	for(int i = 0; i < params.get_vol4d(); ++i) {
		hmc_complex tmp = { (cl_double) i, (cl_double) (params.get_vol4d() - i) };
		h_complex[i] = tmp;
	}
	d_readComplex = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, params.get_vol4d() * sizeof(hmc_complex), h_complex, &err);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	d_float2 = clCreateBuffer(context, CL_MEM_WRITE_ONLY, params.get_vol4d() * sizeof(cl_float2), 0, &err);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	d_complex = clCreateBuffer(context, CL_MEM_WRITE_ONLY, params.get_vol4d() * sizeof(hmc_complex), 0, &err );
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
	cl_int err = clEnqueueReadBuffer(*queue, d_complex, CL_TRUE, 0, params.get_vol4d() * sizeof(hmc_complex), h_complex, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);

	for(int i = 0; i < params.get_vol4d(); ++i) {
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
	cl_float2 * h_float2 = new cl_float2[params.get_vol4d()];
	cl_int err = clEnqueueReadBuffer(*queue, d_float2, CL_TRUE, 0, params.get_vol4d() * sizeof(cl_float2), h_float2, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);

	for(int i = 0; i < params.get_vol4d(); ++i) {
		cl_float2 ref = {(cl_float) i, (cl_float) (params.get_vol4d() - i)};
		verify(ref, h_float2[i]);
	}

	delete[] h_float2;
}

void Dummyfield::runReadKernel()
{
	static_cast<Device*>(opencl_modules[0])->runReadKernel(d_float2, d_readComplex);
}

