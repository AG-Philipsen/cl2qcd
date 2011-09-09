#include "../opencl.h"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SU2 SU3 Extend
#include <boost/test/unit_test.hpp>

#define NUM_ELEMENTS 1024
#define LOCAL_SIZE 128

Random rnd(15);
extern std::string const version;
std::string const version = "0.1";

class Device : public Opencl {

private:
	inputparameters params;
	cl_kernel extendKernel;
	cl_mem in, out;
	Matrixsu2 * h_in;
	Matrixsu3 * h_out;
	cl_int * h_rand;
	cl_mem d_rand;

	void verify(hmc_complex, hmc_complex);

public:
	Device(cl_device_type device_type) : Opencl() {
		Opencl::init(device_type, &params, 0); /* init in body for proper this-pointer */
	};
	virtual void fill_buffers();
	virtual void fill_kernels();
	virtual void clear_buffers();
	virtual void clear_kernels();
	~Device() {
		finalize();
	};

	void runExtendKernel();
	void verify();
};

BOOST_AUTO_TEST_CASE( CPU )
{
	Device dev(CL_DEVICE_TYPE_CPU);
	dev.runExtendKernel();
	dev.verify();
}

BOOST_AUTO_TEST_CASE( GPU )
{
	Device dev(CL_DEVICE_TYPE_GPU);
	dev.runExtendKernel();
	dev.verify();
}

void Device::fill_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	cl_int err;

	h_in = new Matrixsu2[NUM_ELEMENTS];
	for(int i = 0; i < NUM_ELEMENTS; ++i) {
		h_in[i].e00 = hmc_complex_zero;
		h_in[i].e01 = hmc_complex_zero;
		h_in[i].e10 = hmc_complex_zero;
		h_in[i].e11 = hmc_complex_zero;
	}
	BOOST_REQUIRE(h_in);

	in = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, NUM_ELEMENTS * sizeof(Matrixsu2), h_in, &err );
	BOOST_REQUIRE_EQUAL(err,CL_SUCCESS);

	// for simplicity initialize input with 0

	h_out = new Matrixsu3[NUM_ELEMENTS];
	BOOST_REQUIRE(h_out);
	out = clCreateBuffer(context, CL_MEM_WRITE_ONLY, NUM_ELEMENTS * sizeof(Matrixsu3), 0, &err );
	BOOST_REQUIRE_EQUAL(err,CL_SUCCESS);

	h_rand = new cl_int[NUM_ELEMENTS];
	BOOST_REQUIRE(h_rand);
	for(int i = 0; i < NUM_ELEMENTS; ++i) {
		h_rand[i] = (i%3)+1; // high quality random numbers ;)
	}
	d_rand = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, NUM_ELEMENTS * sizeof(cl_int), h_rand, &err );
	BOOST_REQUIRE_EQUAL(err,CL_SUCCESS);

	return;
}

void Device::fill_kernels()
{
	Opencl::fill_kernels();

	extendKernel = createKernel("extendKernel") << basic_opencl_code << "tests/su2su3extend.cl";
}

void Device::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	clReleaseMemObject(in);
	clReleaseMemObject(out);

	delete[] h_in;
	delete[] h_out;

	return;
}

void Device::clear_kernels()
{
	Opencl::clear_kernels();

	clReleaseKernel(extendKernel);

	return;
}

void Device::runExtendKernel() {
	cl_int err;
	cl_ulong elems = NUM_ELEMENTS;
	err = clSetKernelArg(extendKernel, 0, sizeof(cl_mem), &out);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(extendKernel, 1, sizeof(cl_mem), &in);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(extendKernel, 2, sizeof(cl_mem), &d_rand);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(extendKernel, 3, sizeof(cl_ulong), &elems);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	enqueueKernel(extendKernel, NUM_ELEMENTS, LOCAL_SIZE);
}

void Device::verify(hmc_complex left, hmc_complex right) {
	BOOST_REQUIRE_EQUAL(left.re, right.re);
	BOOST_REQUIRE_EQUAL(left.im, right.im);
}


void Device::verify() {
	// get stuff from device
	cl_int err = clEnqueueReadBuffer(queue, out, CL_TRUE, 0, NUM_ELEMENTS * sizeof(Matrixsu3), h_out, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);

	for(size_t i = 0; i < NUM_ELEMENTS; ++i) {
		const int rand = h_rand[i];
		const Matrixsu3 m = h_out[i];

		verify(m.e00, (rand==2) ? hmc_complex_one : hmc_complex_zero);
		verify(m.e01, hmc_complex_zero);
		verify(m.e02, hmc_complex_zero);
		verify(m.e10, hmc_complex_zero);
		verify(m.e11, (rand==3) ? hmc_complex_one : hmc_complex_zero);
		verify(m.e12, hmc_complex_zero);
#ifndef _RECONSTRUCT_TWELVE_
		verify(m.e20, hmc_complex_zero);
		verify(m.e21, hmc_complex_zero);
		verify(m.e22, (rand==1) ? hmc_complex_one : hmc_complex_zero);
#endif
	}
}
