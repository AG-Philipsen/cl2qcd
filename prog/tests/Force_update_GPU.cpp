#include "../opencl_module_hmc.h"
#include "../gaugefield_hybrid.h"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE force_update
#include <boost/test/unit_test.hpp>

Random rnd(15);
extern std::string const version;
std::string const version = "0.1";

class Device : public Opencl_Module_Hmc {

	cl_kernel testKernel;
public:
	Device(cl_command_queue queue, inputparameters* params, int maxcomp, string double_ext) : Opencl_Module_Hmc() {
		Opencl_Module_Hmc::init(queue, 0, params, maxcomp, double_ext); /* init in body for proper this-pointer */
	};
	~Device() {
		finalize();
	};

	void runTestKernel(cl_mem in, cl_mem out, int gs, int ls);
	void fill_kernels();
	void clear_kernels();
};

class Dummyfield : public Gaugefield_hybrid {

public:
	Dummyfield(cl_device_type device_type) : Gaugefield_hybrid() {
		std::stringstream tmp;
#ifdef _USEDOUBLEPREC_
		tmp << SOURCEDIR << "/tests/f_gauge_input_1";
#else
		tmp << SOURCEDIR << "/tests/f_gauge_input_1_single";
#endif
		params.readfile(tmp.str().c_str());

		init(1, device_type, &params);
	};

	virtual void init_tasks();
	virtual void finalize_opencl();

	hmc_float get_squarenorm(int which);
	void verify(hmc_float, hmc_float);
	void verify_result(hmc_float, hmc_float);
	void runTestKernel();

private:
	void fill_buffers();
	void clear_buffers();
	inputparameters params;
	cl_mem in, out;
	cl_mem sqnorm;
	hmc_float * gm_in;
	hmc_float * gm_out;
	Matrixsu3 * gf_in;

};

BOOST_AUTO_TEST_CASE( F_UPDATE )
{
	logger.info() << "Init CPU device";
	//params.print_info_inverter("m_gpu");
	Dummyfield cpu(CL_DEVICE_TYPE_CPU);
	logger.info() << "gaugeobservables: ";
	cpu.print_gaugeobservables_from_task(0, 0);
	logger.info() << "|in|^2:";
	hmc_float cpu_back = cpu.get_squarenorm(0);
	logger.info() << "|out|^2:";
	hmc_float cpu_back2 = cpu.get_squarenorm(1);
	cpu.runTestKernel();
	logger.info() << "|out - 0.12*in|^2:";
	hmc_float cpu_res;
	cpu_res = cpu.get_squarenorm(1);
	BOOST_MESSAGE("Tested CPU");

	logger.info() << "Init GPU device";
	//params.print_info_inverter("m_gpu");
	Dummyfield dummy(CL_DEVICE_TYPE_GPU);
	logger.info() << "gaugeobservables: ";
	dummy.print_gaugeobservables_from_task(0, 0);
	logger.info() << "|in|^2:";
	hmc_float gpu_back = dummy.get_squarenorm(0);
	logger.info() << "|out|^2:";
	hmc_float gpu_back2 = dummy.get_squarenorm(1);
	dummy.runTestKernel();
	logger.info() << "|out + 0.12*in |^2:";
	hmc_float gpu_res;
	gpu_res = dummy.get_squarenorm(1);
	BOOST_MESSAGE("Tested GPU");

	logger.info() << "Compare CPU and GPU results";
	logger.info() << "Input vectors:";
	cpu.verify(cpu_back, gpu_back);
	cpu.verify(cpu_back2, gpu_back2);
	logger.info() << "Output vectors:";
	cpu.verify_result(cpu_res, gpu_res);
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

void fill_with_one(hmc_float * sf_in, int size)
{
	for(int i = 0; i < size; ++i) {
		sf_in[i] = 1.;
	}
	return;
}

void fill_with_random(hmc_float * sf_in, int size, int switcher)
{
	if(switcher == 1) {
		Random rnd_loc(123456);
		for(int i = 0; i < size; ++i) {
			sf_in[i] = rnd_loc.doub();
		}
	} else if (switcher == 2) {
		Random rnd_loc(789101);
		for(int i = 0; i < size; ++i) {
			sf_in[i] = rnd_loc.doub();
		}
	}
	return;
}


void Dummyfield::fill_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	cl_int err;

	cl_context context = opencl_modules[0]->get_context();

	int NUM_ELEMENTS_AE = params.get_gaugemomentasize() * params.get_su3algebrasize();

	gm_in = new hmc_float[NUM_ELEMENTS_AE];
	gm_out = new hmc_float[NUM_ELEMENTS_AE];

	//use the variable use_cg to switch between cold and random input sf
	if(get_parameters()->get_use_cg() == true) {
		fill_with_one(gm_in, NUM_ELEMENTS_AE);
		fill_with_one(gm_out, NUM_ELEMENTS_AE);
	} else {
		fill_with_random(gm_in, NUM_ELEMENTS_AE, 1);
		fill_with_random(gm_out, NUM_ELEMENTS_AE, 2);
	}
	BOOST_REQUIRE(gm_in);
	BOOST_REQUIRE(gm_out);

	Device * device = static_cast<Device*>(opencl_modules[0]);

	size_t ae_buf_size = device->get_gaugemomentum_buffer_size();
	//create buffer for sf on device (and copy sf_in to both for convenience)

	in = clCreateBuffer(context, CL_MEM_READ_WRITE , ae_buf_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	device->importGaugemomentumBuffer(in, reinterpret_cast<ae*>(gm_in));
	out = clCreateBuffer(context, CL_MEM_READ_WRITE, ae_buf_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	device->importGaugemomentumBuffer(out, reinterpret_cast<ae*>(gm_out));

	//create buffer for squarenorm on device
	sqnorm = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float), 0, &err);
}

void Device::fill_kernels()
{
	//one only needs some kernels up to now. to save time during compiling they are put in here by hand
	Opencl_Module_Hmc::fill_kernels();

	testKernel = createKernel("md_update_gaugemomenta") << basic_fermion_code << "types_hmc.h"  << "operations_gaugemomentum.cl" << "md_update_gaugemomenta.cl";
}

void Dummyfield::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	clReleaseMemObject(in);
	clReleaseMemObject(out);
	clReleaseMemObject(sqnorm);

	delete[] gm_in;
	delete[] gm_out;
}

void Device::clear_kernels()
{
	clReleaseKernel(testKernel);
	Opencl_Module::clear_kernels();
}

void Device::runTestKernel(cl_mem out, cl_mem in, int gs, int ls)
{
	cl_int err;
	hmc_float eps = 0.12;
	err = clSetKernelArg(testKernel, 0, sizeof(hmc_float), &eps);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 1, sizeof(cl_mem), &out);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 2, sizeof(cl_mem), &in);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);

	enqueueKernel(testKernel, gs, ls);
}

hmc_float Dummyfield::get_squarenorm(int which)
{
	//which controlls if the in or out-vector is looked at
	if(which == 0) static_cast<Device*>(opencl_modules[0])->set_float_to_gaugemomentum_squarenorm_device(in, sqnorm);
	if(which == 1) static_cast<Device*>(opencl_modules[0])->set_float_to_gaugemomentum_squarenorm_device(out, sqnorm);
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
		logger.info() << "cpu: " << cpu << "\tgpu: " << gpu;
		BOOST_REQUIRE_EQUAL(1, 0);
	}

}

void Dummyfield::verify_result(hmc_float cpu, hmc_float gpu)
{
	//this is too much required, since rounding errors can occur
	//  BOOST_REQUIRE_EQUAL(cpu, gpu);
	//instead, test if the two number agree within some percent
	hmc_float dev = (cpu - gpu) / cpu / 100.;
	if(dev < 1e-10)
		logger.info() << "CPU and GPU result agree within accuary of " << 1e-10;
	else {
		logger.info() << "CPU and GPU result DO NOT agree within accuary of " << 1e-10;
		logger.info() << "cpu: " << cpu << "\tgpu: " << gpu;
		BOOST_REQUIRE_EQUAL(1, 0);
	}
	//in addition, calculate the result again on the host
	hmc_float host_res = 0.;
	for(int i = 0; i < params.get_gaugemomentasize()*params.get_su3algebrasize(); i++) {
		hmc_float tmp = gm_out[i] + gm_in[i] * 0.12;
		host_res += tmp * tmp;
	}
	logger.info() << "the actual result computed on the host is:\t" << host_res;
	dev = (cpu - host_res) / cpu / 100.;
	if(dev < 1e-10)
		logger.info() << "CPU and GPU result agree with the actual result within accuary of " << 1e-10;
	else {
		logger.info() << "CPU and GPU result DO NOT agree with the actual result within accuary of " << 1e-10;
		logger.info() << "cpu: " << cpu << "\tgpu: " << gpu;
		BOOST_REQUIRE_EQUAL(1, 0);
	}


}

void Dummyfield::runTestKernel()
{
	int gs = 0, ls = 0;
	if(opencl_modules[0]->get_device_type() == CL_DEVICE_TYPE_GPU) {
		gs = get_parameters()->get_spinorfieldsize();
		ls = 64;
	} else if(opencl_modules[0]->get_device_type() == CL_DEVICE_TYPE_CPU) {
		gs = opencl_modules[0]->get_max_compute_units();
		ls = 1;
	}
	static_cast<Device*>(opencl_modules[0])->runTestKernel(out, in, gs, ls);
}

