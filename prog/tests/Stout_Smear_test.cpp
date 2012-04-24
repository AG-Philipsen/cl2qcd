#include "../opencl_module.h"
#include "../gaugefield_hybrid.h"

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

//global parameters
hmc_float rho = 0.01;
int iter = 1;


class Device : public Opencl_Module {

	cl_kernel testKernel;
public:
	Device(cl_command_queue queue, inputparameters* params, int maxcomp, string double_ext) : Opencl_Module() {
		Opencl_Module::init(queue, params, maxcomp, double_ext); /* init in body for proper this-pointer */
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
	Dummyfield(cl_device_type device_type) : Gaugefield_hybrid() {
		std::stringstream tmp;
#ifdef _USEDOUBLEPREC_
		tmp << SOURCEDIR << "/tests/stout_smear_test_double";
#else
		tmp << SOURCEDIR << "/tests/stout_smear_test_single";
#endif
		params.readfile(tmp.str().c_str());

		init(1, device_type, &params);
	};

	virtual void init_tasks();
	virtual void finalize_opencl();

	void get_gaugeobservables_from_task(int ntask, hmc_float * plaq, hmc_float * tplaq, hmc_float * splaq, hmc_complex * pol);
	void get_gaugeobservables_from_task(int dummy, int ntask, hmc_float * plaq, hmc_float * tplaq, hmc_float * splaq, hmc_complex * pol);
	void runTestKernel();

private:
	void fill_buffers();
	void clear_buffers();
	inputparameters params;
	cl_mem out;
};

BOOST_AUTO_TEST_CASE( STAPLE_TEST )
{

	hmc_float plaq_cpu, tplaq_cpu, splaq_cpu;
	hmc_complex pol_cpu;
	hmc_float plaq_gpu, tplaq_gpu, splaq_gpu;
	hmc_complex pol_gpu;

	logger.info() << "Init CPU device";
	Dummyfield dummy(CL_DEVICE_TYPE_CPU);
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
	BOOST_MESSAGE("Tested CPU");

	logger.info() << "Init GPU device";
	Dummyfield gpu(CL_DEVICE_TYPE_GPU);
	logger.info() << "gaugeobservables of in field before: ";
	gpu.print_gaugeobservables_from_task(0, 0);
	logger.info() << "gaugeobservables of out field before: ";
	gpu.get_gaugeobservables_from_task(0, 0, &plaq_gpu, &tplaq_gpu, &splaq_gpu, &pol_gpu);
	logger.info() << "plaq: " << plaq_gpu << "\t" << tplaq_gpu  << "\t" << splaq_gpu  << "\t" << pol_gpu.re  << "\t" << pol_gpu.im ;
	gpu.runTestKernel();
	logger.info() << "gaugeobservables of in field after: ";
	gpu.print_gaugeobservables_from_task(0, 0);
	logger.info() << "gaugeobservables of out field after: ";
	gpu.get_gaugeobservables_from_task(0, 0, &plaq_gpu, &tplaq_gpu, &splaq_gpu, &pol_gpu);
	logger.info() << "plaq: " << plaq_gpu << "\t" << tplaq_gpu  << "\t" << splaq_gpu  << "\t" << pol_gpu.re  << "\t" << pol_gpu.im ;
	BOOST_MESSAGE("Tested GPU");

	logger.info() << "test results:";
	logger.info() << "plaq: CPU: " << plaq_cpu << "\tGPU: " << plaq_gpu;
	BOOST_MESSAGE(plaq_cpu << ' ' << plaq_gpu);
	BOOST_CHECK_CLOSE(plaq_cpu, plaq_gpu, 1e-8);
	BOOST_CHECK_CLOSE(tplaq_cpu, tplaq_gpu, 1e-8);
	BOOST_CHECK_CLOSE(splaq_cpu, splaq_gpu, 1e-8);
	CLX_CHECK_CLOSE(pol_cpu, pol_gpu, 1e-8);

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

	int NUM_ELEMENTS = get_num_gaugefield_elems();//params.get_vol4d() * NDIM;


	Matrixsu3 *  gf_tmp  = new Matrixsu3[get_num_gaugefield_elems()];
	//fill tmp gf with ones
	set_gaugefield_cold(gf_tmp);

	Device * device = static_cast<Device*>(opencl_modules[0]);

	out = device->create_rw_buffer(device->getGaugefieldBufferSize());

	//copy cold tmp gf to the device
	device->importGaugefield(out, gf_tmp);

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
	testKernel = createKernel("stout_smear") << basic_opencl_code  <<  "tests/operations_gaugemomentum.cl" << "stout_smear.cl";
}

void Dummyfield::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers
	clReleaseMemObject(out);
}

void Device::clear_kernels()
{
	clReleaseKernel(testKernel);
	Opencl_Module::clear_kernels();
}

void Device::runTestKernel(cl_mem gf, cl_mem out, int gs, int ls)
{
	cl_int err;
	err = clSetKernelArg(testKernel, 0, sizeof(int), &iter  );
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 1, sizeof(cl_mem), &gf);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 2, sizeof(cl_mem), &out);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);

	enqueueKernel(testKernel, gs, ls);
	clFinish(get_queue() );
}

void Dummyfield::runTestKernel()
{
	int gs = 0, ls = 0;
	if(opencl_modules[0]->get_device_type() == CL_DEVICE_TYPE_GPU) {
		gs = get_parameters()->get_vol4d();
		ls = 64;
	} else {
		gs = opencl_modules[0]->get_max_compute_units();
		ls = 1;
	}
	Device * device = static_cast<Device*>(opencl_modules[0]);
	device->runTestKernel(device->get_gaugefield(), out, gs, ls);

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
	opencl_modules[ntask]->gaugeobservables(out, plaq, tplaq, splaq, pol);
}
