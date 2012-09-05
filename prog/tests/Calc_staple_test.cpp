#include "../opencl_module_hmc.h"
#include "../gaugefield_hybrid.h"

#include "../meta/util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE staple_test
#include <boost/test/unit_test.hpp>

extern std::string const version;
std::string const version = "0.1";

#define CLX_CHECK_CLOSE(left, right, precision) \
{ \
  BOOST_CHECK_CLOSE(left.re, right.re, precision); \
  BOOST_CHECK_CLOSE(left.im, right.im, precision); \
}

class Device : public Opencl_Module {

	cl_kernel testKernel;
public:
	Device(cl_command_queue queue, const meta::Inputparameters& params, int maxcomp, std::string double_ext, unsigned int dev_rank) : Opencl_Module(params) {
		Opencl_Module::init(queue, maxcomp, double_ext, dev_rank); /* init in body for proper this-pointer */
	};
	~Device() {
		finalize();
	};

	void runTestKernel(cl_mem gf, cl_mem out, int gs, int ls);
	void fill_kernels();
	void clear_kernels();
};

class Dummyfield : public Gaugefield_hybrid {

public:
	Dummyfield(cl_device_type device_type, const meta::Inputparameters& params)
		: Gaugefield_hybrid(params) {
		init(1, device_type);
	};

	virtual void init_tasks();
	virtual void finalize_opencl();

	hmc_float runTestKernel();

private:
	void fill_buffers();
	void clear_buffers();
	cl_mem out;
	hmc_float * host_out;

};

BOOST_AUTO_TEST_CASE( STAPLE_TEST )
{
	logger.info() << "Init CPU device";
	//params.print_info_inverter("m_gpu");
	// reset RNG
	prng_init(13);
	std::string src = std::string(SOURCEDIR) + "/tests/" + "staple_input_1";
	const char* _params_cpu[] = {"foo", src.c_str(), "--use_gpu=false"};
	meta::Inputparameters params_cpu(3, _params_cpu);
	Dummyfield cpu(CL_DEVICE_TYPE_CPU, params_cpu);
	logger.info() << "gaugeobservables: ";
	cpu.print_gaugeobservables_from_task(0, 0);
	hmc_float cpu_back = cpu.runTestKernel();
	BOOST_CHECK_CLOSE(-1.39070784162e+02, cpu_back, 1e-8);
	BOOST_MESSAGE("Tested CPU");

	logger.info() << "Init GPU device";
	//params.print_info_inverter("m_gpu");
	// reset RNG
	prng_init(13);
	src = std::string(SOURCEDIR) + "/tests/" + "staple_input_1";
	const char* _params_gpu[] = {"foo", src.c_str(), "--use_gpu=true"};
	meta::Inputparameters params_gpu = meta::Inputparameters(3, _params_gpu);
	Dummyfield dummy(CL_DEVICE_TYPE_GPU, params_gpu);
	logger.info() << "gaugeobservables: ";
	dummy.print_gaugeobservables_from_task(0, 0);
	hmc_float gpu_back = dummy.runTestKernel();
	BOOST_CHECK_CLOSE(-1.39070784162e+02, gpu_back, 1e-8);
	//u_res = dummy.get_squarenorm(1);
	BOOST_MESSAGE("Tested GPU");

	logger.info() << "cpu: " << std::scientific << std::setprecision(11) << cpu_back << "\tgpu: " << gpu_back;

	BOOST_CHECK_CLOSE(cpu_back, gpu_back, 1e-8);

	//CP: in case of a cold config, the result is calculable easily
	src = std::string(SOURCEDIR) + "/tests/" + "staple_input_cold1";
	const char* _params_cpu_cold[] = {"foo", src.c_str(), "--use_gpu=false"};
	meta::Inputparameters params_cpu_cold = meta::Inputparameters(3, _params_cpu_cold);
	Dummyfield cpu_cold(CL_DEVICE_TYPE_CPU, params_cpu_cold);
	hmc_float cold_back = cpu_cold.runTestKernel();
	hmc_float cold_ref = 18 * NDIM * meta::get_vol4d(dummy.get_parameters());
	BOOST_CHECK_CLOSE(cold_back, cold_ref, 1e-8);

	// TODO test further input files, especially larger sizes and anisotropic case
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

	int NUM_ELEMENTS = meta::get_vol4d(get_parameters());

	host_out = new hmc_float[NUM_ELEMENTS];
	BOOST_REQUIRE(host_out);

	size_t buf_size = NUM_ELEMENTS * sizeof(hmc_float);
	out = clCreateBuffer(context, CL_MEM_READ_ONLY , buf_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

}

void Device::fill_kernels()
{
	//one only needs some kernels up to now. to save time during compiling they are put in here by hand
	Opencl_Module::fill_kernels();

	//to this end, one has to set the needed files by hand
	basic_opencl_code = ClSourcePackage() << "opencl_header.cl" << "operations_geometry.cl" << "operations_complex.cl"
	                    << "operations_matrix_su3.cl" << "operations_matrix.cl" << "operations_gaugefield.cl";

	testKernel = createKernel("staple_test") << basic_opencl_code  << "/tests/staple_test.cl";

}

void Dummyfield::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	clReleaseMemObject(out);

	delete[] host_out;
}

void Device::clear_kernels()
{
	clReleaseKernel(testKernel);
	Opencl_Module::clear_kernels();
}

void Device::runTestKernel(cl_mem gf, cl_mem out, int gs, int ls)
{
	cl_int err;
	err = clSetKernelArg(testKernel, 0, sizeof(cl_mem), &gf);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 1, sizeof(cl_mem), &out );
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);

	enqueueKernel(testKernel, gs, ls);
}

hmc_float Dummyfield::runTestKernel()
{
	hmc_float res = 0;
	int gs, ls;
	if(opencl_modules[0]->get_device_type() == CL_DEVICE_TYPE_GPU) {
		gs = meta::get_vol4d(get_parameters());
		ls = 64;
	} else {
		gs = opencl_modules[0]->get_max_compute_units();
		ls = 1;
	}
	Device * device = static_cast<Device*>(opencl_modules[0]);
	device->runTestKernel(device->get_gaugefield(), out, gs, ls);

	int NUM_ELEMENTS = meta::get_vol4d(get_parameters());
	//copy the result of the kernel to host
	size_t size = NUM_ELEMENTS * sizeof(hmc_float);
	cl_int clerr = clEnqueueReadBuffer(opencl_modules[0]->get_queue(), out, CL_TRUE, 0, size, host_out, 0, NULL, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clEnqueueReadBuffer", __FILE__, __LINE__);

	//sum up all elements in the result buffer
	for(int i = 0; i < NUM_ELEMENTS; i++) {
		res += host_out[i];
	}
	return res;
}

