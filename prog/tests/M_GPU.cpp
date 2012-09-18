#include "../opencl_module_fermions.h"
#include "../gaugefield_hybrid.h"
#include "../meta/util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Fermionmatrix
#include <boost/test/unit_test.hpp>

extern std::string const version;
std::string const version = "0.1";

class Device : public Opencl_Module_Fermions {

	cl_kernel testKernel;

public:
	Device(cl_command_queue queue, const meta::Inputparameters& params, int maxcomp, std::string double_ext, unsigned int dev_rank) : Opencl_Module_Fermions(params) {
		Opencl_Module_Fermions::init(queue, maxcomp, double_ext, dev_rank); /* init in body for proper this-pointer */
	};
	~Device() {
		finalize();
	};

	void runTestKernel(cl_mem in, cl_mem out, int gs, int ls, hmc_float kappa, hmc_float mubar);
	void fill_kernels();
	void clear_kernels();
};

const std::string SOURCEFILE = std::string(SOURCEDIR) + "/tests/m_gpu_input_1";
const char * PARAMS[] = {"foo", SOURCEFILE.c_str()};
const meta::Inputparameters INPUT(2, PARAMS);

class Dummyfield : public Gaugefield_hybrid {

public:
	Dummyfield(cl_device_type device_type, const hardware::System * system) : Gaugefield_hybrid(system) {
		init(1, device_type);
	};

	virtual void init_tasks();
	virtual void finalize_opencl();

	hmc_float get_squarenorm(int which);
	void verify(hmc_float, hmc_float);
	void runTestKernel();

private:
	void fill_buffers();
	void clear_buffers();
	cl_mem in, out;
	cl_mem sqnorm;
	spinor * sf_in;
	spinor * sf_out;
};

BOOST_AUTO_TEST_CASE( M_TM )
{
	logger.info() << "Init CPU device";
	//params.print_info_inverter("m_gpu");
	hardware::System system_cpu(INPUT);
	Dummyfield cpu(CL_DEVICE_TYPE_CPU, &system_cpu);
	logger.info() << "gaugeobservables: ";
	cpu.print_gaugeobservables_from_task(0, 0);
	logger.info() << "|phi|^2:";
	hmc_float cpu_back = cpu.get_squarenorm(0);
	cpu.runTestKernel();
	logger.info() << "|M phi|^2:";
	hmc_float cpu_res;
	cpu.print_gaugeobservables_from_task(0, 0);
	BOOST_MESSAGE("Tested CPU");

	logger.info() << "Init GPU device";
	//params.print_info_inverter("m_gpu");
	hardware::System system_gpu(INPUT);
	Dummyfield dummy(CL_DEVICE_TYPE_GPU, &system_gpu);
	logger.info() << "gaugeobservables: ";
	dummy.print_gaugeobservables_from_task(0, 0);
	logger.info() << "|phi|^2:";
	hmc_float gpu_back = dummy.get_squarenorm(0);
	dummy.runTestKernel();
	logger.info() << "|M phi|^2:";
	hmc_float gpu_res;
	gpu_res = dummy.get_squarenorm(1);
	BOOST_MESSAGE("Tested GPU");

	logger.info() << "Compare CPU and GPU results";
	logger.info() << "Input vectors:";
	cpu.verify(cpu_back, gpu_back);
	logger.info() << "Output vectors:";
	cpu.verify(cpu_res, gpu_res);
}

// BOOST_AUTO_TEST_CASE( CPU )
// {
//  logger.info() << "Init dummy device";
//  //params.print_info_inverter("m_gpu");
//  Dummyfield dummy(CL_DEVICE_TYPE_CPU);
//  logger.info() << "gaugeobservables: ";
//  dummy.print_gaugeobservables_from_task(0, 0);
//  logger.info() << "|phi|^2:";
//  dummy.verify(0);
//  dummy.runTestKernel();
//  logger.info() << "|M phi|^2:";
//  dummy.verify(1);
//  BOOST_MESSAGE("Tested CPU");
// }
//
// BOOST_AUTO_TEST_CASE( GPU )
// {
//  logger.info() << "Init dummy device";
//  //params.print_info_inverter("m_gpu");
//  Dummyfield dummy(CL_DEVICE_TYPE_GPU);
//  logger.info() << "gaugeobservables: ";
//  dummy.print_gaugeobservables_from_task(0, 0);
//  logger.info() << "|phi|^2:";
//  dummy.verify(0);
//  dummy.runTestKernel();
//  logger.info() << "|M phi|^2:";
//  dummy.verify(1);
//  BOOST_MESSAGE("Tested GPU");
// }

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

void fill_sf_with_one(spinor * sf_in, int size)
{
	for(int i = 0; i < size; ++i) {
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
	return;
}

void fill_sf_with_random(spinor * sf_in, int size)
{
	prng_init(123456);
	for(int i = 0; i < size; ++i) {
		sf_in[i].e0.e0.re = prng_double();
		sf_in[i].e0.e1.re = prng_double();
		sf_in[i].e0.e2.re = prng_double();
		sf_in[i].e1.e0.re = prng_double();
		sf_in[i].e1.e1.re = prng_double();
		sf_in[i].e1.e2.re = prng_double();
		sf_in[i].e2.e0.re = prng_double();
		sf_in[i].e2.e1.re = prng_double();
		sf_in[i].e2.e2.re = prng_double();
		sf_in[i].e3.e0.re = prng_double();
		sf_in[i].e3.e1.re = prng_double();
		sf_in[i].e3.e2.re = prng_double();

		sf_in[i].e0.e0.im = prng_double();
		sf_in[i].e0.e1.im = prng_double();
		sf_in[i].e0.e2.im = prng_double();
		sf_in[i].e1.e0.im = prng_double();
		sf_in[i].e1.e1.im = prng_double();
		sf_in[i].e1.e2.im = prng_double();
		sf_in[i].e2.e0.im = prng_double();
		sf_in[i].e2.e1.im = prng_double();
		sf_in[i].e2.e2.im = prng_double();
		sf_in[i].e3.e0.im = prng_double();
		sf_in[i].e3.e1.im = prng_double();
		sf_in[i].e3.e2.im = prng_double();
	}
	return;
}


void Dummyfield::fill_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	cl_int err;

	cl_context context = opencl_modules[0]->get_context();

	size_t NUM_ELEMENTS_SF =  meta::get_spinorfieldsize(get_parameters());

	sf_in = new spinor[NUM_ELEMENTS_SF];
	sf_out = new spinor[NUM_ELEMENTS_SF];

	//use the variable use_cg to switch between cold and random input sf
	if(get_parameters().get_solver() == meta::Inputparameters::cg) fill_sf_with_one(sf_in, NUM_ELEMENTS_SF);
	else fill_sf_with_random(sf_in, NUM_ELEMENTS_SF);
	BOOST_REQUIRE(sf_in);

	size_t sf_buf_size = NUM_ELEMENTS_SF * sizeof(spinor);
	logger.warn() << sf_buf_size;
	//create buffer for sf on device (and copy sf_in to both for convenience)
	in = clCreateBuffer(context, CL_MEM_READ_WRITE, sf_buf_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	err = clEnqueueWriteBuffer(static_cast<Device*>(opencl_modules[0])->get_queue(), in, CL_TRUE, 0, sf_buf_size, sf_in, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	out = clCreateBuffer(context, CL_MEM_READ_WRITE, sf_buf_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	err = clEnqueueWriteBuffer(static_cast<Device*>(opencl_modules[0])->get_queue(), out, CL_TRUE, 0, sf_buf_size, sf_in, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	//create buffer for squarenorm on device
	sqnorm = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float), 0, &err);
}

void Device::fill_kernels()
{
	Opencl_Module_Fermions::fill_kernels();

	testKernel = createKernel("M_tm_plus") << basic_fermion_code  << "fermionmatrix.cl" << "fermionmatrix_m_tm_plus.cl";

}

void Dummyfield::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	clReleaseMemObject(in);
	clReleaseMemObject(out);
	clReleaseMemObject(sqnorm);

	delete[] sf_in;
	delete[] sf_out;
}

void Device::clear_kernels()
{
	clReleaseKernel(testKernel);
	Opencl_Module::clear_kernels();
}

void Device::runTestKernel(cl_mem out, cl_mem in, int gs, int ls, hmc_float kappa, hmc_float mubar)
{
	cl_mem gf = get_gaugefield();
	cl_int err;
	err = clSetKernelArg(testKernel, 0, sizeof(cl_mem), &in);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 1, sizeof(cl_mem), &gf);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 2, sizeof(cl_mem), &out);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 3, sizeof(hmc_float), &kappa);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 4, sizeof(hmc_float), &mubar);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);

	enqueueKernel(testKernel, gs, ls);
}

hmc_float Dummyfield::get_squarenorm(int which)
{
	//which controlls if the in or out-vector is looked at
	if(which == 0) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_device(in, sqnorm);
	if(which == 1) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_device(out, sqnorm);
	// get stuff from device
	hmc_float result;
	cl_int err = clEnqueueReadBuffer(*queue, sqnorm, CL_TRUE, 0, sizeof(hmc_float), &result, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	logger.info() << result;

	return result;
}

void Dummyfield::verify(hmc_float cpu, hmc_float gpu)
{
	//this is too much required, since rounding errors can occur
	//  BOOST_REQUIRE_EQUAL(cpu, gpu);
	//instead, test if the two number agree within some percent
	hmc_float dev = (cpu - gpu) / cpu / 100.;
	if(dev < 1e-10)
		logger.info() << "CPU and GPU result agree within accuary of " << 1e-10;
	else {
		logger.info() << "CPU and GPU result DO NOT agree within accuary of " << 1e-10;
		BOOST_REQUIRE_EQUAL(cpu, gpu);
	}
}

void Dummyfield::runTestKernel()
{
	int gs = 0, ls = 0;
	if(opencl_modules[0]->get_device_type() == CL_DEVICE_TYPE_GPU) {
		gs = meta::get_spinorfieldsize(get_parameters());
		ls = 64;
	} else if(opencl_modules[0]->get_device_type() == CL_DEVICE_TYPE_CPU) {
		gs = opencl_modules[0]->get_max_compute_units();
		ls = 1;
	}
	Device * device = static_cast<Device*>(opencl_modules[0]);
	device->runTestKernel(out, in, gs, ls, get_parameters().get_kappa(), meta::get_mubar(get_parameters()));
}

