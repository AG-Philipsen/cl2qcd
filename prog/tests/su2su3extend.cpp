#include "../opencl_module.h"
#include "../gaugefield_hybrid.h"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SU2 SU3 Extend
#include <boost/test/unit_test.hpp>

#define NUM_ELEMENTS 1024
#define LOCAL_SIZE 128

extern std::string const version;
std::string const version = "0.1";

class Device : public Opencl_Module {

	cl_kernel extendKernel;

public:
	Device(cl_command_queue queue, inputparameters* params, int maxcomp, string double_ext, unsigned int dev_rank) : Opencl_Module() {
		Opencl_Module::init(queue, params, maxcomp, double_ext, dev_rank); /* init in body for proper this-pointer */
	};
	~Device() {
		finalize();
	};

	void runExtendKernel(cl_mem out, cl_mem in, cl_mem d_rand, cl_ulong elems);
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

	void verify();
	void runExtendKernel();

private:
	void verify(hmc_complex, hmc_complex);
	void fill_buffers();
	void clear_buffers();
	inputparameters params;
	Matrixsu2 * h_in;
	Matrixsu3 * h_out;
	cl_int * h_rand;
	cl_mem in, out;
	cl_mem d_rand;

};

BOOST_AUTO_TEST_CASE( CPU )
{
	Dummyfield dummy(CL_DEVICE_TYPE_CPU);
	dummy.runExtendKernel();
	dummy.verify();
	BOOST_MESSAGE("Tested CPU");
}

BOOST_AUTO_TEST_CASE( GPU )
{
	Dummyfield dummy(CL_DEVICE_TYPE_GPU);
	dummy.runExtendKernel();
	dummy.verify();
	BOOST_MESSAGE("Tested GPU");
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
	// don't invoke parent function as we don't require the original buffers

	cl_int err;

	cl_context context = opencl_modules[0]->get_context();

	h_in = new Matrixsu2[NUM_ELEMENTS];
	for(int i = 0; i < NUM_ELEMENTS; ++i) {
		h_in[i].e00 = hmc_complex_zero;
		h_in[i].e01 = hmc_complex_zero;
		h_in[i].e10 = hmc_complex_zero;
		h_in[i].e11 = hmc_complex_zero;
	}
	BOOST_REQUIRE(h_in);

	in = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, NUM_ELEMENTS * sizeof(Matrixsu2), h_in, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	// for simplicity initialize input with 0

	h_out = new Matrixsu3[NUM_ELEMENTS];
	BOOST_REQUIRE(h_out);
	out = clCreateBuffer(context, CL_MEM_WRITE_ONLY, NUM_ELEMENTS * sizeof(Matrixsu3), 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	h_rand = new cl_int[NUM_ELEMENTS];
	BOOST_REQUIRE(h_rand);
	for(int i = 0; i < NUM_ELEMENTS; ++i) {
		h_rand[i] = (i % 3) + 1; // high quality random numbers ;)
	}
	d_rand = clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, NUM_ELEMENTS * sizeof(cl_int), h_rand, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
}

void Device::fill_kernels()
{
	Opencl_Module::fill_kernels();

	extendKernel = createKernel("extendKernel") << basic_opencl_code << "tests/su2su3extend.cl";
}

void Dummyfield::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	clReleaseMemObject(in);
	clReleaseMemObject(out);

	delete[] h_in;
	delete[] h_out;
}

void Device::clear_kernels()
{
	clReleaseKernel(extendKernel);
	Opencl_Module::clear_kernels();
}

void Device::runExtendKernel(cl_mem out, cl_mem in, cl_mem d_rand, cl_ulong elems)
{
	cl_int err;
	err = clSetKernelArg(extendKernel, 0, sizeof(cl_mem), &out);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(extendKernel, 1, sizeof(cl_mem), &in);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(extendKernel, 2, sizeof(cl_mem), &d_rand);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(extendKernel, 3, sizeof(cl_ulong), &elems);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	enqueueKernel(extendKernel, NUM_ELEMENTS);
}

void Dummyfield::verify(hmc_complex left, hmc_complex right)
{
	BOOST_REQUIRE_EQUAL(left.re, right.re);
	BOOST_REQUIRE_EQUAL(left.im, right.im);
}


void Dummyfield::verify()
{
	// get stuff from device
	cl_int err = clEnqueueReadBuffer(*queue, out, CL_TRUE, 0, NUM_ELEMENTS * sizeof(Matrixsu3), h_out, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);

	for(size_t i = 0; i < NUM_ELEMENTS; ++i) {
		const int rand = h_rand[i];
		const Matrixsu3 m = h_out[i];

		verify(m.e00, (rand == 2) ? hmc_complex_one : hmc_complex_zero);
		verify(m.e01, hmc_complex_zero);
		verify(m.e02, hmc_complex_zero);
		verify(m.e10, hmc_complex_zero);
		verify(m.e11, (rand == 3) ? hmc_complex_one : hmc_complex_zero);
		verify(m.e12, hmc_complex_zero);
		verify(m.e20, hmc_complex_zero);
		verify(m.e21, hmc_complex_zero);
		verify(m.e22, (rand == 1) ? hmc_complex_one : hmc_complex_zero);
	}
}

void Dummyfield::runExtendKernel()
{
	static_cast<Device*>(opencl_modules[0])->runExtendKernel(out, in, d_rand, NUM_ELEMENTS);
}
