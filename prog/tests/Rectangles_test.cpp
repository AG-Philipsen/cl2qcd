#include "../opencl_module_hmc.h"
#include "../gaugefield_hybrid.h"
#include "../meta/util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Rectangles
#include <boost/test/unit_test.hpp>

extern std::string const version;
std::string const version = "0.1";
std::string const exec_name = "f_tlsym";

class Device : public Opencl_Module_Hmc {

	//for rectangles
	cl_kernel testKernel;
	cl_kernel testKernel_b;
	//for force_gauge_tlsym
	cl_kernel testKernel2;
	meta::Counter counter1, counter2, counter3, counter4;
public:
	Device(cl_command_queue queue, const meta::Inputparameters& params, int maxcomp, std::string double_ext, unsigned int dev_rank) : Opencl_Module_Hmc(params, &counter1, &counter2, &counter3, &counter4) {
		Opencl_Module_Hmc::init(queue, maxcomp, double_ext, dev_rank); /* init in body for proper this-pointer */
	};
	~Device() {
		finalize();
	};

	//for rectangles
	void runTestKernel(cl_mem out, cl_mem gf, int gs, int ls);
	//for force_gauge_tlsym
	void runTestKernel2(cl_mem out, cl_mem gf, int gs, int ls);
	void fill_kernels();
	void clear_kernels();
};

const std::string SOURCEFILE = std::string(SOURCEDIR)
#ifdef _USEDOUBLEPREC_
                               + "/tests/f_rectangles_input_1";
#else
                               + "/tests/f_rectangles_input_1_single";
#endif
const char * PARAMS[] = {"foo", SOURCEFILE.c_str()};
const meta::Inputparameters INPUT(2, PARAMS);

class Dummyfield : public Gaugefield_hybrid {

public:
  Dummyfield(cl_device_type device_type, meta::Inputparameters inputfile) : Gaugefield_hybrid(inputfile) {
    init(1, device_type);
    meta::print_info_hmc(exec_name.c_str(), inputfile);
  };
  
  virtual void init_tasks();
  virtual void finalize_opencl();

	hmc_float get_squarenorm();
	hmc_float get_rect();
  bool verify(hmc_float, hmc_float, hmc_float);
	void runTestKernel();
	void runTestKernel2();

private:
	void fill_buffers();
	void clear_buffers();
	cl_mem out;
	cl_mem sqnorm;
	cl_mem rect_value;
	Matrixsu3 * gf_in;
	hmc_float * sf_out;

};

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

void fill_with_zero(hmc_float * sf_in, int size)
{
	for(int i = 0; i < size; ++i) {
		sf_in[i] = 0.;
	}
	return;
}

void Dummyfield::fill_buffers()
{
	// don't invoke parent function as we don't require the original buffers
	cl_int err;
	cl_context context = opencl_modules[0]->get_context();
	int NUM_ELEMENTS_AE = meta::get_vol4d(get_parameters()) * NDIM * meta::get_su3algebrasize();
	sf_out = new hmc_float[NUM_ELEMENTS_AE];
	fill_with_zero(sf_out, NUM_ELEMENTS_AE);

	Device * device = static_cast<Device*>(opencl_modules[0]);
	out = clCreateBuffer(context, CL_MEM_WRITE_ONLY, device->get_gaugemomentum_buffer_size(), 0, &err);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	device->importGaugemomentumBuffer(out, reinterpret_cast<ae*>(sf_out));

	//create buffer for squarenorm on device
	sqnorm = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float), 0, &err);
	//create buffer for value of rectangles
	rect_value = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float), 0, &err);
}

void Device::fill_kernels()
{
	//one only needs some kernels up to now. to save time during compiling they are put in here by hand
	Opencl_Module_Hmc::fill_kernels();

	testKernel = createKernel("rectangles") << basic_fermion_code << "gaugeobservables_rectangles.cl";
	testKernel_b = createKernel("rectangles_reduction") << basic_opencl_code << "gaugeobservables_rectangles.cl";
	//at the time of writing this kernel, the OpenCL compiler crashed the kernel using optimizations
	testKernel2 = createKernel("gauge_force_tlsym") << basic_fermion_code << "types_hmc.h"  << "operations_gaugemomentum.cl" << "force_gauge_tlsym.cl";
}

void Dummyfield::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers
	clReleaseMemObject(out);
	clReleaseMemObject(sqnorm);
	clReleaseMemObject(rect_value);
	delete[] sf_out;
}

void Device::clear_kernels()
{
	clReleaseKernel(testKernel);
	clReleaseKernel(testKernel2);
	Opencl_Module::clear_kernels();
}

void Device::runTestKernel(cl_mem out, cl_mem gf, int gs, int ls)
{
	logger.trace() << "start with testKernel";
	cl_int err;
	cl_uint num_groups = (gs + ls - 1) / ls;
	int global_buf_size_float = sizeof(hmc_float) * num_groups;
	cl_mem clmem_rect_buf_glob = 0;
	if( clmem_rect_buf_glob == 0 ) clmem_rect_buf_glob = create_rw_buffer(global_buf_size_float);
	int buf_loc_size_float = sizeof(hmc_float) * ls;

	//set arguments
	// run local rectangles calculation and first part of reduction
	int clerr = clSetKernelArg(testKernel, 0, sizeof(cl_mem), &gf);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(testKernel, 1, sizeof(cl_mem), &clmem_rect_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(testKernel, 2, buf_loc_size_float, NULL);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(testKernel, gs, ls);

	// run second part of rectangles reduction

	clerr = clSetKernelArg(testKernel_b, 0, sizeof(cl_mem), &clmem_rect_buf_glob);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(testKernel_b, 1, sizeof(cl_mem), &out);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);
	clerr = clSetKernelArg(testKernel_b, 2, sizeof(cl_uint), &num_groups);
	if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

	enqueueKernel(testKernel_b, 1, 1);
	clReleaseMemObject(clmem_rect_buf_glob);
	logger.trace() << "done with testKernel";

}

void Device::runTestKernel2(cl_mem out, cl_mem gf, int gs, int ls)
{
	logger.trace() << "start with testKernel2";
	cl_int err;
	err = clSetKernelArg(testKernel2, 0, sizeof(cl_mem), &gf);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel2, 1, sizeof(cl_mem), &out);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	enqueueKernel(testKernel2, gs, ls);
	logger.trace() << "done with testKernel2";
}


hmc_float Dummyfield::get_squarenorm()
{
	static_cast<Device*>(opencl_modules[0])->set_float_to_gaugemomentum_squarenorm_device(out, sqnorm);
	// get stuff from device
	hmc_float result;
	cl_int err = clEnqueueReadBuffer(*queue, sqnorm, CL_TRUE, 0, sizeof(hmc_float), &result, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	return result;
}

hmc_float Dummyfield::get_rect()
{
	// get stuff from device
	hmc_float result;
	cl_int err = clEnqueueReadBuffer(*queue, rect_value, CL_TRUE, 0, sizeof(hmc_float), &result, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	//logger.info() << result;
	//hmc_float norm =  NDIM * (NDIM - 1) * NC *  meta::get_vol4d(get_parameters());
	//logger.info() << "in the correct normalization (12*VOL*NC = " << norm << ") this reads:";
	//logger.info() << result / norm;
	return result;
}

bool Dummyfield::verify(hmc_float cpu, hmc_float gpu, hmc_float prec)
{
	//this is too much required, since rounding errors can occur
	//  BOOST_REQUIRE_EQUAL(cpu, gpu);
	//instead, test if the two number agree up to some precision prec
	hmc_float dev = (abs(cpu) - abs(gpu)) / cpu;
	if(abs(dev) < prec) {
		return true;
	} else {
		return false;
	}
}

void Dummyfield::runTestKernel()
{
	int gs = 0, ls = 0;
	if(opencl_modules[0]->get_device_type() == CL_DEVICE_TYPE_GPU) {
		gs = get_spinorfieldsize(get_parameters());
		ls = 64;
	} else if(opencl_modules[0]->get_device_type() == CL_DEVICE_TYPE_CPU) {
		gs = opencl_modules[0]->get_max_compute_units();
		ls = 1;
	}
	Device * device = static_cast<Device*>(opencl_modules[0]);
	device->runTestKernel(rect_value, device->get_gaugefield(), gs, ls);
}

void Dummyfield::runTestKernel2()
{
	int gs = 0, ls = 0;
	if(opencl_modules[0]->get_device_type() == CL_DEVICE_TYPE_GPU) {
		gs = meta::get_spinorfieldsize(get_parameters());
		ls = 128;
	} else if(opencl_modules[0]->get_device_type() == CL_DEVICE_TYPE_CPU) {
		gs = opencl_modules[0]->get_max_compute_units();
		ls = 1;
	}
	Device * device = static_cast<Device*>(opencl_modules[0]);
	device->runTestKernel2(out, device->get_gaugefield(), gs, ls);
}


BOOST_AUTO_TEST_CASE( F_GAUGE )
{
  logger.info() << "Test CPU and GPU version of kernel";
  logger.info() << "\tf_tlsym";
  logger.info() << "and";
  logger.info() << "\trectangles";
  logger.info() << "against reference value";
  logger.info() << "(this is done in one test right now because both kernels use similar functions...";
  
  logger.info() << "Init CPU device";
  Dummyfield cpu(CL_DEVICE_TYPE_CPU, INPUT);
  logger.info() << "gaugeobservables: ";
  cpu.print_gaugeobservables_from_task(0, 0);
  logger.info() << "calc rectangles value:";
  cpu.runTestKernel();
  hmc_float cpu_rect = cpu.get_rect();
  logger.info() << cpu_rect;
  cpu.runTestKernel2();
  logger.info() << "|f_gauge|^2:";
  hmc_float cpu_res;
  cpu_res = cpu.get_squarenorm();
  logger.info() << cpu_res;
  BOOST_MESSAGE("Tested CPU");

  logger.info() << "Init GPU device";
  Dummyfield gpu(CL_DEVICE_TYPE_GPU, INPUT);
  logger.info() << "gaugeobservables: ";
  gpu.print_gaugeobservables_from_task(0, 0);
  logger.info() << "calc rectangles value:";
  gpu.runTestKernel();
  hmc_float gpu_rect = gpu.get_rect();
  logger.info() << gpu_rect;
  gpu.runTestKernel2();
  logger.info() << "|f_gauge|^2:";
  hmc_float gpu_res;
  gpu_res = gpu.get_squarenorm();
  logger.info() << gpu_res;
  BOOST_MESSAGE("Tested GPU");

  logger.info() << "Choosing reference value";
  //this is the value of the force measured in tmlqcd for the same config
  hmc_float force_ref_value = 4080.6745694080;
  hmc_float rect_ref_value = 1103.2398401620;
  logger.info() << "reference value force:\t" << force_ref_value;
  logger.info() << "reference value rect:\t" << rect_ref_value;

  hmc_float prec = 1e-10;  
  logger.info() << "acceptance precision: " << prec;
  logger.info() << "Compare CPU result to reference value";
  logger.info() << "rect:";
  bool res1a = cpu.verify(cpu_rect, rect_ref_value, prec);
  if(res1a) {
    logger.info() << "CPU and reference value agree within accuary of " << prec;
  } else {
    logger.info() << "CPU and reference value DO NOT agree within accuary of " << prec;
    BOOST_REQUIRE_EQUAL(1, 0);
  }
  logger.info() << "force:";
  bool res1b = cpu.verify(cpu_res, force_ref_value, prec);
  if(res1b) {
    logger.info() << "CPU and reference value agree within accuary of " << prec;
  } else {
    logger.info() << "CPU and reference value DO NOT agree within accuary of " << prec;
    BOOST_REQUIRE_EQUAL(1, 0);
  }

  logger.info() << "Compare GPU result to reference value";
  logger.info() << "rect:";
  bool res2a = gpu.verify(gpu_rect, rect_ref_value, prec);
  if(res2a) {
    logger.info() << "GPU and reference value agree within accuary of " << prec;
  } else {
    logger.info() << "GPU and reference value DO NOT agree within accuary of " << prec;
    BOOST_REQUIRE_EQUAL(1, 0);
  }
  logger.info() << "force:";
  bool res2b = gpu.verify(gpu_res, rect_ref_value, prec);
  if(res2b) {
    logger.info() << "GPU and reference value agree within accuary of " << prec;
  } else {
    logger.info() << "GPU and reference value DO NOT agree within accuary of " << prec;
    BOOST_REQUIRE_EQUAL(1, 0);
  }

  logger.info() << "Compare CPU and GPU results";
  logger.info() << "rect:";
  bool res3a = cpu.verify(cpu_rect, gpu_rect, prec);
  if(res3a) {
    logger.info() << "CPU and GPU result agree within accuary of " << prec;
  } else {
    logger.info() << "CPU and GPU result DO NOT agree within accuary of " << prec;
    BOOST_REQUIRE_EQUAL(1, 0);
  }
  logger.info() << "force:";
  bool res3b = cpu.verify(cpu_res, gpu_res, prec);
  if(res3b) {
    logger.info() << "CPU and GPU result agree within accuary of " << prec;
  } else {
    logger.info() << "CPU and GPU result DO NOT agree within accuary of " << prec;
    BOOST_REQUIRE_EQUAL(1, 0);
  }
}

