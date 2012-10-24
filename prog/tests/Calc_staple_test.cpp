#include "../opencl_module_hmc.h"
#include "../gaugefield_hybrid.h"

#include "../meta/util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE staple_test
#include <boost/test/unit_test.hpp>

extern std::string const version;
std::string const version = "0.1";
std::string const exec_name = "staple_test";

class Device : public Opencl_Module_Gaugefield {

	cl_kernel testKernel;
public:
	Device(const meta::Inputparameters& params, hardware::Device * device) : Opencl_Module_Gaugefield(params, device) {
		fill_kernels();
	};
	~Device() {
		clear_kernels();
	};

	void runTestKernel(const hardware::buffers::SU3 * gf, const hardware::buffers::Plain<hmc_float> * out, int gs, int ls);
	void fill_kernels();
	void clear_kernels();
};


class Dummyfield : public Gaugefield_hybrid {
public:
	Dummyfield(const hardware::System * system) : Gaugefield_hybrid(system) {
		auto inputfile = system->get_inputparameters();
		init(1, inputfile.get_use_gpu() ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU);
		meta::print_info_hmc(exec_name.c_str(), inputfile);
	};
	virtual void init_tasks();
	virtual void finalize_opencl();

	hmc_float get_squarenorm();
	hmc_float runTestKernel();

private:
	void fill_buffers();
	void clear_buffers();
	const hardware::buffers::Plain<hmc_float> * out;
	hmc_float * host_out;
};

void Dummyfield::init_tasks()
{
	opencl_modules = new Opencl_Module* [get_num_tasks()];
	opencl_modules[0] = new Device(get_parameters(), get_device_for_task(0));

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
	int NUM_ELEMENTS = meta::get_vol4d(get_parameters());
	host_out = new hmc_float[NUM_ELEMENTS];
	BOOST_REQUIRE(host_out);

	out = new hardware::buffers::Plain<hmc_float>(NUM_ELEMENTS, opencl_modules[0]->get_device());
}

void Device::fill_kernels()
{
	testKernel = createKernel("staple_test") << get_device()->get_gaugefield_code()->get_sources() << "/tests/staple_test.cl";
}

void Dummyfield::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers
	delete out;
	delete[] host_out;
}

void Device::clear_kernels()
{
	clReleaseKernel(testKernel);
}

void Device::runTestKernel(const hardware::buffers::SU3 * gf, const hardware::buffers::Plain<hmc_float> * out, int gs, int ls)
{
	cl_int err;
	err = clSetKernelArg(testKernel, 0, sizeof(cl_mem), gf->get_cl_buffer());
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 1, sizeof(cl_mem), out->get_cl_buffer());
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);

	get_device()->enqueue_kernel(testKernel, gs, ls);
}

hmc_float Dummyfield::runTestKernel()
{
	hmc_float res = 0;
	int gs, ls;
	if(get_device_for_task(0)->get_device_type() == CL_DEVICE_TYPE_GPU) {
		gs = meta::get_vol4d(get_parameters());
		ls = 64;
	} else {
		gs = get_device_for_task(0)->get_num_compute_units();
		ls = 1;
	}
	Device * device = static_cast<Device*>(opencl_modules[0]);
	device->runTestKernel(device->get_gaugefield(), out, gs, ls);

	int NUM_ELEMENTS = meta::get_vol4d(get_parameters());
	//copy the result of the kernel to host
	out->dump(host_out);

	//sum up all elements in the result buffer
	for(int i = 0; i < NUM_ELEMENTS; i++) {
		res += host_out[i];
	}
	return res;
}


BOOST_AUTO_TEST_CASE( STAPLE_TEST )
{
	logger.info() << "Test kernel";
	logger.info() << "\tcalc_staple";
	logger.info() << "against reference value";

	int param_expect = 4;
	logger.info() << "expect parameters:";
	logger.info() << "\texec_name\tinputfile\tgpu_usage\trec12_usage";
	//get number of parameters
	int num_par = boost::unit_test::framework::master_test_suite().argc;
	if(num_par < param_expect) {
		logger.fatal() << "need more inputparameters! Got only " << num_par << ", expected " << param_expect << "! Aborting...";
		exit(-1);
	}

	//get input file that has been passed as an argument
	const char*  inputfile =  boost::unit_test::framework::master_test_suite().argv[1];
	logger.info() << "inputfile used: " << inputfile;
	//get use_gpu = true/false that has been passed as an argument
	const char*  gpu_opt =  boost::unit_test::framework::master_test_suite().argv[2];
	logger.info() << "GPU usage: " << gpu_opt;
	//get use_rec12 = true/false that has been passed as an argument
	const char* rec12_opt =  boost::unit_test::framework::master_test_suite().argv[3];
	logger.info() << "rec12 usage: " << rec12_opt;

	logger.info() << "Init device";
	const char* _params_cpu[] = {"foo", inputfile, gpu_opt, rec12_opt};
	meta::Inputparameters params(param_expect, _params_cpu);
	hardware::System system(params);
	Dummyfield cpu(&system);
	logger.info() << "gaugeobservables: ";
	cpu.print_gaugeobservables_from_task(0, 0);
	logger.info() << "Run kernel";
	logger.info() << "running test kernel";
	hmc_float cpu_res = cpu.runTestKernel();
	logger.info() << "result:";
	logger.info() << cpu_res;

	logger.info() << "Choosing reference value and acceptance precision";
	hmc_float ref_val = params.get_test_ref_value();
	logger.info() << "reference value:\t" << ref_val;
	hmc_float prec = params.get_solver_prec();
	logger.info() << "acceptance precision: " << prec;

	logger.info() << "Compare result to reference value";
	BOOST_REQUIRE_CLOSE(cpu_res, ref_val, prec);
	logger.info() << "Done";
	BOOST_MESSAGE("Test done");
}

