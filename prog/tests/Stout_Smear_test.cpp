#include "../opencl_module.h"
#include "../gaugefield_hybrid.h"
#include "../meta/util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE stout_smear_test
#include <boost/test/unit_test.hpp>

extern std::string const version;
std::string const version = "0.1";
std::string const exec_name = "stout_smear_test";

#define CLX_CHECK_CLOSE(left, right, precision) \
{ \
  BOOST_CHECK_CLOSE(left.re, right.re, precision); \
  BOOST_CHECK_CLOSE(left.im, right.im, precision); \
}

//global parameters
hmc_float rho = 0.01;
int iter = 1;

class Device : public Opencl_Module {
	cl_kernel testKernel;
public:
	Device(const meta::Inputparameters& params, hardware::Device * device) : Opencl_Module(params, device), out(get_gaugefield()->get_elements(), device) {
		Opencl_Module::init(); /* init in body for proper this-pointer */
	};
	~Device() {
		finalize();
	};

	void runTestKernel(const hardware::buffers::SU3 * gf, const hardware::buffers::SU3 * out, int gs, int ls);
	void fill_kernels();
	void clear_kernels();

	const hardware::buffers::SU3 out;
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

	void get_gaugeobservables_from_task(int ntask, hmc_float * plaq, hmc_float * tplaq, hmc_float * splaq, hmc_complex * pol);
	void get_gaugeobservables_from_task(int dummy, int ntask, hmc_float * plaq, hmc_float * tplaq, hmc_float * splaq, hmc_complex * pol);
	void runTestKernel();

private:
	void fill_buffers();
};

void Dummyfield::init_tasks()
{
	opencl_modules = new Opencl_Module* [get_num_tasks()];
	opencl_modules[0] = new Device(get_parameters(), get_device_for_task(0));

	fill_buffers();
}

void Dummyfield::finalize_opencl()
{
	Gaugefield_hybrid::finalize_opencl();
}

void Dummyfield::fill_buffers()
{
	Matrixsu3 *  gf_tmp  = new Matrixsu3[get_num_gaugefield_elems()];
	//fill tmp gf with ones
	set_gaugefield_cold(gf_tmp);

	Device * device = static_cast<Device*>(opencl_modules[0]);

	//copy cold tmp gf to the device
	device->importGaugefield(&device->out, gf_tmp);

	delete[] gf_tmp;

}

void Device::fill_kernels()
{
	//one only needs some kernels up to now. to save time during compiling they are put in here by hand
	Opencl_Module::fill_kernels();

	//to this end, one has to set the needed files by hand
	basic_opencl_code = ClSourcePackage() << "opencl_header.cl" << "operations_geometry.cl" << "operations_complex.cl"
	                    << "operations_matrix_su3.cl" << "operations_matrix.cl" << "operations_gaugefield.cl";

	//this has to be the same as in opencl_module
	// in fact, the kernel has already been build in the above call!!
	//  testKernel = createKernel("stout_smear") << basic_opencl_code  <<  "tests/operations_gaugemomentum.cl" << "stout_smear.cl";
}
void Device::clear_kernels()
{
	clReleaseKernel(testKernel);
	Opencl_Module::clear_kernels();
}

void Device::runTestKernel(const hardware::buffers::SU3 * gf, const hardware::buffers::SU3 * out, int gs, int ls)
{
	cl_int err;
	err = clSetKernelArg(testKernel, 0, sizeof(cl_mem), gf->get_cl_buffer());
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 1, sizeof(cl_mem), out->get_cl_buffer());
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);

	get_device()->enqueue_kernel(testKernel, gs, ls);
	clFinish(get_queue() );
}

void Dummyfield::runTestKernel()
{
	//CP: this currently causes a segfault!!!
	Device * device = static_cast<Device*>(opencl_modules[0]);
	static_cast<Device*>(opencl_modules[0])->stout_smear_device( device->get_gaugefield()  , &device->out);

	int gs = 0, ls = 0;
	if(get_device_for_task(0)->get_device_type() == CL_DEVICE_TYPE_GPU) {
		gs = meta::get_vol4d(get_parameters());
		ls = 64;
	} else {
		gs = get_device_for_task(0)->get_num_compute_units();
		ls = 1;
	}
	//  Device * device = static_cast<Device*>(opencl_modules[0]);
	//  device->runTestKernel(device->get_gaugefield(), out, gs, ls);

	return;
}

void Dummyfield::get_gaugeobservables_from_task(int ntask, hmc_float * plaq, hmc_float * tplaq, hmc_float * splaq, hmc_complex * pol)
{
	if( ntask < 0 || ntask > get_num_tasks() ) throw Print_Error_Message("devicetypes index out of range", __FILE__, __LINE__);
	opencl_modules[ntask]->gaugeobservables(plaq, tplaq, splaq, pol);
}

//this is just out of laziness, a copy of the function above
void Dummyfield::get_gaugeobservables_from_task(int dummy, int ntask, hmc_float * plaq, hmc_float * tplaq, hmc_float * splaq, hmc_complex * pol)
{
	dummy = 0;
	if( ntask < 0 || ntask > get_num_tasks() ) throw Print_Error_Message("devicetypes index out of range", __FILE__, __LINE__);
	opencl_modules[ntask]->gaugeobservables(&static_cast<Device*>(opencl_modules[0])->out, plaq, tplaq, splaq, pol);
}

BOOST_AUTO_TEST_CASE( STOUT_SMEAR )
{
	logger.info() << "Test kernel";
	logger.info() << "\tstout_smear";
	logger.info() << "against reference value";

	logger.fatal() << "A segfault appears when the kernel is called using the proper module fct! Exit..";
	BOOST_REQUIRE_EQUAL(1., 0.);

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
	Dummyfield dummy(&system);

	hmc_float plaq_cpu, tplaq_cpu, splaq_cpu;
	hmc_complex pol_cpu;
	hmc_float plaq_gpu, tplaq_gpu, splaq_gpu;
	hmc_complex pol_gpu;

	logger.info() << "gaugeobservables of in field before: ";
	dummy.print_gaugeobservables_from_task(0, 0);
	logger.info() << "gaugeobservables of out field before: ";
	dummy.get_gaugeobservables_from_task(0, 0, &plaq_cpu, &tplaq_cpu, &splaq_cpu, &pol_cpu);
	logger.info() << "plaq: " << plaq_cpu << "\t" << tplaq_cpu  << "\t" << splaq_cpu  << "\t" << pol_cpu.re  << "\t" << pol_cpu.im ;
	dummy.runTestKernel();
	logger.info() << "gaugeobservables of in field after: ";
	dummy.print_gaugeobservables_from_task(0, 0);
	logger.info() << "gaugeobservables of out field after: ";
	dummy.get_gaugeobservables_from_task(0, 0, &plaq_cpu, &tplaq_cpu, &splaq_cpu, &pol_cpu);
	logger.info() << "plaq: " << plaq_cpu << "\t" << tplaq_cpu  << "\t" << splaq_cpu  << "\t" << pol_cpu.re  << "\t" << pol_cpu.im ;

	logger.info() << "Choosing reference value and acceptance precision";
	hmc_float ref_val = params.get_test_ref_value();
	logger.info() << "reference value:\t" << ref_val;
	hmc_float prec = params.get_solver_prec();
	logger.info() << "acceptance precision: " << prec;

	logger.info() << "Compare result to reference value";
	BOOST_REQUIRE_CLOSE(plaq_cpu, ref_val, prec);
	logger.info() << "Done";
	BOOST_MESSAGE("Test done");
}
