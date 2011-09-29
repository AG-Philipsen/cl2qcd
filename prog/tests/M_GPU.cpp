#include "../opencl_module_spinors.h"
#include "../gaugefield_hybrid.h"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SU2 SU3 Extend
#include <boost/test/unit_test.hpp>

Random rnd(15);
extern std::string const version;
std::string const version = "0.1";

class Device : public Opencl_Module_Spinors {

	cl_kernel testKernel;

public:
	Device(cl_command_queue queue, inputparameters* params, int maxcomp, string double_ext) : Opencl_Module_Spinors() {
		Opencl_Module_Spinors::init(queue, 0, params, maxcomp, double_ext); /* init in body for proper this-pointer */
	};
	~Device() {
		finalize();
	};

	void runTestKernel(cl_mem in, cl_mem out, cl_mem gf, int gs, int ls);
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

	void verify(int which);
	void runTestKernel();

private:
	void verify(hmc_complex, hmc_complex);
	void fill_buffers();
	void clear_buffers();
	inputparameters params;
	cl_mem in, out, gf;
	cl_mem sqnorm;
	spinor * sf_in;
	spinor * sf_out;
	Matrixsu3 * gf_in;

};

BOOST_AUTO_TEST_CASE( CPU )
{
	BOOST_MESSAGE("Init dummy device");
	Dummyfield dummy(CL_DEVICE_TYPE_CPU);
	BOOST_MESSAGE("|phi|^2:");
	dummy.verify(0);
	dummy.runTestKernel();
	BOOST_MESSAGE("|M phi|^2:");
	dummy.verify(1);
	BOOST_MESSAGE("Tested CPU");
}

BOOST_AUTO_TEST_CASE( GPU )
{
	BOOST_MESSAGE("Init dummy device");
	Dummyfield dummy(CL_DEVICE_TYPE_GPU);
	BOOST_MESSAGE("|phi|^2:");
	dummy.verify(0);
	dummy.runTestKernel();
	BOOST_MESSAGE("|M phi|^2:");
	dummy.verify(1);
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

	int NUM_ELEMENTS_SF = params.get_spinorfieldsize();
	
	sf_in = new spinor[NUM_ELEMENTS_SF];
	sf_out = new spinor[NUM_ELEMENTS_SF];
	//init spinorfield with 1s
	for(int i = 0; i < NUM_ELEMENTS_SF; ++i) {
		sf_in[i].e0.e0 = hmc_complex_one;
		sf_in[i].e0.e1 = hmc_complex_one;
		sf_in[i].e0.e2 = hmc_complex_one;
		sf_in[i].e1.e0 = hmc_complex_one;
		sf_in[i].e1.e1 = hmc_complex_one;
		sf_in[i].e1.e2 = hmc_complex_one;
		sf_in[i].e2.e0 = hmc_complex_one;
		sf_in[i].e2.e1 = hmc_complex_one;
		sf_in[i].e2.e2 = hmc_complex_one;
		sf_in[i].e3.e0 = hmc_complex_one;
		sf_in[i].e3.e1 = hmc_complex_one;
		sf_in[i].e3.e2 = hmc_complex_one;
	}
	BOOST_REQUIRE(sf_in);

	//create buffer for sf on device (and copy sf_in)
	in = clCreateBuffer(context, CL_MEM_READ_ONLY, params.get_sf_buf_size(), sf_in, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	out = clCreateBuffer(context, CL_MEM_WRITE_ONLY, params.get_sf_buf_size(), 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	
	int NUM_ELEMENTS_GF = params.get_gaugefieldsize();
	gf_in = new Matrixsu3[NUM_ELEMENTS_GF];
	set_gaugefield_cold(gf_in);
	BOOST_REQUIRE(gf_in);
	
	//create buffer for gf on device (and copy gf)
	gf = clCreateBuffer(context, CL_MEM_READ_ONLY, params.get_gf_buf_size(), gf_in, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	
	//create buffer for squarenorm on device
	sqnorm = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float), 0, &err);
}

void Device::fill_kernels()
{
	Opencl_Module_Spinors::fill_kernels();
	testKernel = createKernel("M") << basic_opencl_code << "fermionmatrx_GPU.cl";
}

void Dummyfield::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	clReleaseMemObject(in);
	clReleaseMemObject(out);
	clReleaseMemObject(gf);
	clReleaseMemObject(sqnorm);

	delete[] sf_in;
	delete[] sf_out;
	delete[] gf_in;
}

void Device::clear_kernels()
{
	clReleaseKernel(testKernel);
	Opencl_Module::clear_kernels();
}

void Device::runTestKernel(cl_mem out, cl_mem in, cl_mem gf, int gs, int ls)
{
	cl_int err;
	err = clSetKernelArg(testKernel, 0, sizeof(cl_mem), &in);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 1, sizeof(cl_mem), &out);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 2, sizeof(cl_mem), &gf);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	
	enqueueKernel(testKernel, gs, ls);
}

void Dummyfield::verify(int which)
{
	//which controlls if the in or out-vector is looked at
	if(which == 0) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_device(in, sqnorm);
	if(which == 1) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_device(out, sqnorm);
	// get stuff from device
	hmc_float result;
	cl_int err = clEnqueueReadBuffer(*queue, sqnorm, CL_TRUE, 0, sizeof(hmc_float), &result, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	logger.info() << result;
}

void Dummyfield::runTestKernel()
{
	static_cast<Device*>(opencl_modules[0])->runTestKernel(out, in, gf, params.get_spinorfieldsize(), 1);
}
