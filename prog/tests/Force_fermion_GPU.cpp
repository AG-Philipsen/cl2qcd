#include "../opencl_module_hmc.h"
#include "../gaugefield_hybrid.h"

#include "../meta/util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Fermionforce
#include <boost/test/unit_test.hpp>

extern std::string const version;
std::string const version = "0.1";

class Device : public Opencl_Module_Hmc {

	cl_kernel testKernel;
	meta::Counter counter1, counter2, counter3, counter4;
public:
	Device(cl_command_queue queue, const meta::Inputparameters& params, hardware::Device * device, int maxcomp, std::string double_ext, unsigned int dev_rank) : Opencl_Module_Hmc(params, device, &counter1, &counter2, &counter3, &counter4) {
		Opencl_Module_Hmc::init(queue, maxcomp, double_ext, dev_rank); /* init in body for proper this-pointer */
	};
	~Device() {
		finalize();
	};

	void runTestKernel(cl_mem in1, cl_mem in2, cl_mem out, cl_mem gf, int gs, int ls, hmc_float);
	void fill_kernels();
	void clear_kernels();
};

const std::string SOURCEFILE = std::string(SOURCEDIR)
#ifdef _USEDOUBLEPREC_
                               + "/tests/f_gauge_input_1";
#else
                               + "/tests/f_gauge_input_1_single";
#endif
const char * PARAMS[] = {"foo", SOURCEFILE.c_str()};
const meta::Inputparameters INPUT(2, PARAMS);

class Dummyfield : public Gaugefield_hybrid {

public:
	Dummyfield(cl_device_type device_type, hardware::System * system) : Gaugefield_hybrid(system) {
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
	cl_mem in1, in2, out;
	cl_mem sqnorm;
	spinor * sf_in1;
	spinor * sf_in2;
	ae * ae_out;
	Matrixsu3 * gf_in;

};

BOOST_AUTO_TEST_CASE( F_FERMION )
{
	logger.info() << "Init CPU device";
	//params.print_info_inverter("m_gpu");
	// reset RNG
	prng_init(13);
	hardware::System system_cpu(INPUT);
	Dummyfield cpu(CL_DEVICE_TYPE_CPU, &system_cpu);
	logger.info() << "gaugeobservables: ";
	cpu.print_gaugeobservables_from_task(0, 0);
	logger.info() << "|phi_1|^2:";
	hmc_float cpu_back = cpu.get_squarenorm(0);
	logger.info() << "|phi_2|^2:";
	hmc_float cpu_back2 = cpu.get_squarenorm(1);
	cpu.runTestKernel();
	logger.info() << "|M phi|^2:";
	hmc_float cpu_res = cpu.get_squarenorm(2);
	BOOST_MESSAGE("Tested CPU");

	logger.info() << "Init GPU device";
	//params.print_info_inverter("m_gpu");
	// reset RNG
	prng_init(13);
	hardware::System system_gpu(INPUT);
	Dummyfield dummy(CL_DEVICE_TYPE_GPU, &system_gpu);
	logger.info() << "gaugeobservables: ";
	dummy.print_gaugeobservables_from_task(0, 0);
	logger.info() << "|phi_1|^2:";
	hmc_float gpu_back = dummy.get_squarenorm(0);
	logger.info() << "|phi_2|^2:";
	hmc_float gpu_back2 = dummy.get_squarenorm(1);
	dummy.runTestKernel();
	logger.info() << "|M phi|^2:";
	hmc_float gpu_res = dummy.get_squarenorm(2);
	BOOST_MESSAGE("Tested GPU");

	logger.info() << "Compare CPU and GPU results";
	logger.info() << "Input vectors:";
	cpu.verify(cpu_back, gpu_back);
	cpu.verify(cpu_back2, gpu_back2);
	logger.info() << "Output vectors:";
	cpu.verify(cpu_res, gpu_res);
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

void fill_sf_with_float(spinor * sf_in, int size, hmc_float val)
{
	for(int i = 0; i < size; ++i) {
		sf_in[i].e0.e0.re = val;
		sf_in[i].e0.e1.re = val;
		sf_in[i].e0.e2.re = val;
		sf_in[i].e1.e0.re = val;
		sf_in[i].e1.e1.re = val;
		sf_in[i].e1.e2.re = val;
		sf_in[i].e2.e0.re = val;
		sf_in[i].e2.e1.re = val;
		sf_in[i].e2.e2.re = val;
		sf_in[i].e3.e0.re = val;
		sf_in[i].e3.e1.re = val;
		sf_in[i].e3.e2.re = val;

		sf_in[i].e0.e0.im = val;
		sf_in[i].e0.e1.im = val;
		sf_in[i].e0.e2.im = val;
		sf_in[i].e1.e0.im = val;
		sf_in[i].e1.e1.im = val;
		sf_in[i].e1.e2.im = val;
		sf_in[i].e2.e0.im = val;
		sf_in[i].e2.e1.im = val;
		sf_in[i].e2.e2.im = val;
		sf_in[i].e3.e0.im = val;
		sf_in[i].e3.e1.im = val;
		sf_in[i].e3.e2.im = val;
	}
	return;
}

void fill_sf_with_pos(spinor * sf_in, int size)
{
	hmc_float val;
	for(int i = 0; i < size; ++i) {
		val = i;
		sf_in[i].e0.e0.re = val;
		sf_in[i].e0.e1.re = val;
		sf_in[i].e0.e2.re = val;
		sf_in[i].e1.e0.re = val;
		sf_in[i].e1.e1.re = val;
		sf_in[i].e1.e2.re = val;
		sf_in[i].e2.e0.re = val;
		sf_in[i].e2.e1.re = val;
		sf_in[i].e2.e2.re = val;
		sf_in[i].e3.e0.re = val;
		sf_in[i].e3.e1.re = val;
		sf_in[i].e3.e2.re = val;

		sf_in[i].e0.e0.im = val;
		sf_in[i].e0.e1.im = val;
		sf_in[i].e0.e2.im = val;
		sf_in[i].e1.e0.im = val;
		sf_in[i].e1.e1.im = val;
		sf_in[i].e1.e2.im = val;
		sf_in[i].e2.e0.im = val;
		sf_in[i].e2.e1.im = val;
		sf_in[i].e2.e2.im = val;
		sf_in[i].e3.e0.im = val;
		sf_in[i].e3.e1.im = val;
		sf_in[i].e3.e2.im = val;
	}
	return;
}

void fill_with_one(hmc_float * sf_in, int size)
{
	for(int i = 0; i < size; ++i) {
		sf_in[i] = 1.;
	}
	return;
}

ae make_ae(hmc_float e1, hmc_float e2, hmc_float e3, hmc_float e4,
           hmc_float e5, hmc_float e6, hmc_float e7, hmc_float e8)
{
	ae tmp = {e1, e2, e3, e4, e5, e6, e7, e8};
	return tmp;
};

void fill_with_zero(ae * ae, int size)
{
	for(int i = 0; i < size; ++i) {
		ae[i] = make_ae(0., 0., 0., 0., 0., 0., 0., 0.);
	}
	return;
}



void fill_sf_with_random(spinor * sf_in, int size, int seed)
{
	//  Random rnd_loc(seed);
	for(int i = 0; i < size; ++i) {
		prng_init(seed);
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

	int NUM_ELEMENTS_SF = meta::get_spinorfieldsize(get_parameters());

	int NUM_ELEMENTS_AE = meta::get_vol4d(get_parameters()) * NDIM * meta::get_su3algebrasize();

	sf_in1 = new spinor[NUM_ELEMENTS_SF];
	sf_in2 = new spinor[NUM_ELEMENTS_SF];
	ae_out = new ae[NUM_ELEMENTS_AE];

	//use the variable use_cg to switch between cold and random input sf
	if(get_parameters().get_solver() == meta::Inputparameters::cg) {
		fill_sf_with_one(sf_in1, NUM_ELEMENTS_SF);
		fill_sf_with_one(sf_in2, NUM_ELEMENTS_SF);
	} else {
		fill_sf_with_random(sf_in1, NUM_ELEMENTS_SF, 123456);
		fill_sf_with_random(sf_in2, NUM_ELEMENTS_SF, 789101);
	}

	//fill with zeros
	//hmc_float val = 0.;
	//hmc_float val2 = 0.;
	//fill_sf_with_float(sf_in1, NUM_ELEMENTS_SF, val);
	//fill_sf_with_float(sf_in2, NUM_ELEMENTS_SF, val2);


	//  fill_sf_with_pos(sf_in1,2/* NUM_ELEMENTS_SF*/);
	//fill_sf_with_pos(sf_in2,2/* NUM_ELEMENTS_SF*/);

	BOOST_REQUIRE(sf_in1);
	BOOST_REQUIRE(sf_in2);

	fill_with_zero(ae_out, NUM_ELEMENTS_AE);

	size_t sf_buf_size = meta::get_spinorfieldsize(get_parameters()) * sizeof(spinor);
	//create buffer for sf on device (and copy sf_in to both for convenience)

	Device * device = static_cast<Device*>(opencl_modules[0]);

	in1 = clCreateBuffer(context, CL_MEM_READ_ONLY , sf_buf_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	in2 = clCreateBuffer(context, CL_MEM_READ_ONLY , sf_buf_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	err = clEnqueueWriteBuffer(device->get_queue(), in1, CL_TRUE, 0, sf_buf_size, sf_in1, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	err = clEnqueueWriteBuffer(device->get_queue(), in2, CL_TRUE, 0, sf_buf_size, sf_in2, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	out = clCreateBuffer(context, CL_MEM_READ_WRITE, device->get_gaugemomentum_buffer_size(), 0, &err);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	device->importGaugemomentumBuffer(out, reinterpret_cast<ae*>(ae_out));

	//create buffer for squarenorm on device
	sqnorm = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float), 0, &err);
}

void Device::fill_kernels()
{
	//one only needs some kernels up to now. to save time during compiling they are put in here by hand
	Opencl_Module_Hmc::fill_kernels();

	testKernel = createKernel("fermion_force") << basic_fermion_code << "types_hmc.h"  << "operations_gaugemomentum.cl" << "fermionmatrix.cl" << "force_fermion.cl";

}

void Dummyfield::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	clReleaseMemObject(in1);
	clReleaseMemObject(in2);
	clReleaseMemObject(out);
	clReleaseMemObject(sqnorm);

	delete[] sf_in1;
	delete[] sf_in2;
	delete[] ae_out;
}

void Device::clear_kernels()
{
	clReleaseKernel(testKernel);
	Opencl_Module_Hmc::clear_kernels();
}

void Device::runTestKernel(cl_mem out, cl_mem in1, cl_mem in2, cl_mem gf, int gs, int ls, hmc_float kappa)
{
	cl_int err;
	err = clSetKernelArg(testKernel, 0, sizeof(cl_mem), &gf);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 1, sizeof(cl_mem), &in1);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 2, sizeof(cl_mem), &in2);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 3, sizeof(cl_mem), &out);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 4, sizeof(hmc_float), &kappa);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);

	enqueueKernel(testKernel, gs, ls);
}

hmc_float Dummyfield::get_squarenorm(int which)
{
	//which controlls if the in or out-vector is looked at
	if(which == 0) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_device(in1, sqnorm);
	if(which == 1) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_device(in2, sqnorm);
	if(which == 2) static_cast<Device*>(opencl_modules[0])->set_float_to_gaugemomentum_squarenorm_device(out, sqnorm);
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
	if(abs(dev) < 1e-10) {
		logger.info() << "CPU and GPU result agree within accuary of " << 1e-10;
		logger.info() << "cpu: " << cpu << "\tgpu: " << gpu;
	} else {
		logger.info() << "CPU and GPU result DO NOT agree within accuary of " << 1e-10;
		logger.info() << "cpu: " << cpu << "\tgpu: " << gpu;
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
	device->runTestKernel(out, in1, in2, device->get_gaugefield(), gs, ls, get_parameters().get_kappa());
}
