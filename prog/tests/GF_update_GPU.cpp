#include "../opencl_module_hmc.h"
#include "../gaugefield_hybrid.h"
#include "../meta/util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE gf_update
#include <boost/test/unit_test.hpp>

extern std::string const version;
std::string const version = "0.1";

#define CLX_CHECK_CLOSE(left, right, precision) \
{ \
  BOOST_CHECK_CLOSE(left.re, right.re, precision); \
  BOOST_CHECK_CLOSE(left.im, right.im, precision); \
}

class Device : public Opencl_Module_Hmc {

	cl_kernel testKernel;
	meta::Counter counter1, counter2, counter3, counter4;
public:
	Device(cl_command_queue queue, const meta::Inputparameters& params, int maxcomp, std::string double_ext, unsigned int dev_rank) : Opencl_Module_Hmc(params, &counter1, &counter2, &counter3, &counter4) {
		Opencl_Module_Hmc::init(queue, maxcomp, double_ext, dev_rank); /* init in body for proper this-pointer */
	};
	~Device() {
		finalize();
	};

	void runTestKernel(cl_mem in, cl_mem gf, int gs, int ls);
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
	Dummyfield(cl_device_type device_type, const hardware::System * system) : Gaugefield_hybrid(system) {
		init(1, device_type);
	};

	virtual void init_tasks();
	virtual void finalize_opencl();

	hmc_float get_squarenorm(int which);
	void get_gaugeobservables_from_task(int ntask, hmc_float * plaq, hmc_float * tplaq, hmc_float * splaq, hmc_complex * pol);
	void runTestKernel();

private:
	void fill_buffers();
	void clear_buffers();
	cl_mem in, out;
	cl_mem sqnorm;
	hmc_float * gm_in;
	Matrixsu3 * gf_in;

};

BOOST_AUTO_TEST_CASE( F_UPDATE )
{
	logger.info() << "Init CPU device";
	//params.print_info_inverter("m_gpu");
	// reset RNG
	hardware::System system_cpu(INPUT);
	prng_init(13);
	Dummyfield cpu(CL_DEVICE_TYPE_CPU, &system_cpu);
	logger.info() << "gaugeobservables: ";
	cpu.print_gaugeobservables_from_task(0, 0);
	logger.info() << "|in|^2:";
	hmc_float cpu_back = cpu.get_squarenorm(0);
	cpu.runTestKernel();
	logger.info() << "gaugeobservables: ";
	cpu.print_gaugeobservables_from_task(0, 0);

	hmc_float plaq_cpu, tplaq_cpu, splaq_cpu;
	hmc_complex pol_cpu;
	cpu.get_gaugeobservables_from_task(0, &plaq_cpu, &tplaq_cpu, &splaq_cpu, &pol_cpu);

	//u_res = cpu.get_squarenorm(1);
	BOOST_MESSAGE("Tested CPU");

	logger.info() << "Init GPU device";
	//params.print_info_inverter("m_gpu");
	// reset RNG
	prng_init(13);
	hardware::System system_gpu(INPUT);
	Dummyfield dummy(CL_DEVICE_TYPE_GPU, &system_gpu);
	logger.info() << "gaugeobservables: ";
	dummy.print_gaugeobservables_from_task(0, 0);
	logger.info() << "|in|^2:";
	hmc_float gpu_back = dummy.get_squarenorm(0);
	dummy.runTestKernel();
	logger.info() << "gaugeobservables: ";
	dummy.print_gaugeobservables_from_task(0, 0);

	hmc_float plaq_gpu, tplaq_gpu, splaq_gpu;
	hmc_complex pol_gpu;
	dummy.get_gaugeobservables_from_task(0, &plaq_gpu, &tplaq_gpu, &splaq_gpu, &pol_gpu);

	//c_float gpu_res;
	//u_res = dummy.get_squarenorm(1);
	BOOST_MESSAGE("Tested GPU");
	/*
	logger.info() << "Compare CPU and GPU results";
	logger.info() << "Input vectors:";
	cpu.verify(cpu_back, gpu_back);
	cpu.verify(cpu_back2, gpu_back2);
	logger.info() << "Output vectors:";
	cpu.verify(cpu_res, gpu_res);
	*/
	BOOST_MESSAGE(cpu_back << ' ' << gpu_back);
	BOOST_MESSAGE(plaq_cpu << ' ' << plaq_gpu);
	BOOST_CHECK_CLOSE(cpu_back, gpu_back, 1e-8);
	BOOST_CHECK_CLOSE(plaq_cpu, plaq_gpu, 1e-8);
	BOOST_CHECK_CLOSE(tplaq_cpu, tplaq_gpu, 1e-8);
	BOOST_CHECK_CLOSE(splaq_cpu, splaq_gpu, 1e-8);
	CLX_CHECK_CLOSE(pol_cpu, pol_gpu, 1e-8);
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
		prng_init(123456);
		for(int i = 0; i < size; ++i) {
			sf_in[i] = prng_double();
		}
	} else if (switcher == 2) {
		prng_init(789101);
		for(int i = 0; i < size; ++i) {
			sf_in[i] = prng_double();
		}
	}
	return;
}


void Dummyfield::fill_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	cl_int err;

	cl_context context = opencl_modules[0]->get_context();

	int NUM_ELEMENTS_AE = meta::get_vol4d(get_parameters()) * NDIM * meta::get_su3algebrasize();

	gm_in = new hmc_float[NUM_ELEMENTS_AE];

	//use the variable use_cg to switch between cold and random input sf
	if(get_parameters().get_solver() == meta::Inputparameters::cg) {
		fill_with_one(gm_in, NUM_ELEMENTS_AE);
	} else {
		fill_with_random(gm_in, NUM_ELEMENTS_AE, 1);
	}
	BOOST_REQUIRE(gm_in);

	Device * device = static_cast<Device*>(opencl_modules[0]);
	//create buffer for sf on device (and copy sf_in to both for convenience)

	in = clCreateBuffer(context, CL_MEM_READ_ONLY, device->get_gaugemomentum_buffer_size(), 0, &err);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	device->importGaugemomentumBuffer(in, reinterpret_cast<ae*>(gm_in));

	//create buffer for squarenorm on device
	sqnorm = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float), 0, &err);
}

void Device::fill_kernels()
{
	//one only needs some kernels up to now. to save time during compiling they are put in here by hand
	Opencl_Module_Hmc::fill_kernels();

	testKernel = createKernel("md_update_gaugefield") << basic_fermion_code << "types_hmc.h" << "operations_gaugemomentum.cl" << "md_update_gaugefield.cl";
}

void Dummyfield::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	clReleaseMemObject(in);
	clReleaseMemObject(sqnorm);

	delete[] gm_in;
}

void Device::clear_kernels()
{
	clReleaseKernel(testKernel);
	Opencl_Module::clear_kernels();
}

void Device::runTestKernel(cl_mem in, cl_mem gf, int gs, int ls)
{
	cl_int err;
	hmc_float eps = .12;
	err = clSetKernelArg(testKernel, 0, sizeof(hmc_float), &eps);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 1, sizeof(cl_mem), &in );
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 2, sizeof(cl_mem), &gf);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);

	enqueueKernel(testKernel, gs, ls);
}

hmc_float Dummyfield::get_squarenorm(int which)
{
	//which controlls if the in or out-vector is looked at
	if(which == 0) static_cast<Device*>(opencl_modules[0])->set_float_to_gaugemomentum_squarenorm_device(in, sqnorm);
	// get stuff from device
	hmc_float result;
	cl_int err = clEnqueueReadBuffer(*queue, sqnorm, CL_TRUE, 0, sizeof(hmc_float), &result, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	logger.info() << result;
	return result;
}

void Dummyfield::runTestKernel()
{
	int gs, ls;
	if(opencl_modules[0]->get_device_type() == CL_DEVICE_TYPE_GPU) {
		gs = meta::get_spinorfieldsize(get_parameters());
		ls = 64;
	} else {
		gs = opencl_modules[0]->get_max_compute_units();
		ls = 1;
	}
	Device * device = static_cast<Device*>(opencl_modules[0]);
	device->runTestKernel(in, device->get_gaugefield(), gs, ls);
}


void Dummyfield::get_gaugeobservables_from_task(int ntask, hmc_float * plaq, hmc_float * tplaq, hmc_float * splaq, hmc_complex * pol)
{
	if( ntask < 0 || ntask > get_num_tasks() ) throw Print_Error_Message("devicetypes index out of range", __FILE__, __LINE__);
	opencl_modules[ntask]->gaugeobservables(plaq, tplaq, splaq, pol);
}

