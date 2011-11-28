#include "../opencl_module_hmc.h"
#include "../gaugefield_hybrid.h"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Fermionforce_eo
#include <boost/test/unit_test.hpp>

Random rnd(15);
extern std::string const version;
std::string const version = "0.1";

class Device : public Opencl_Module_Hmc {

	cl_kernel testKernel;
  cl_kernel ae_sqn;
public:
	Device(cl_command_queue queue, inputparameters* params, int maxcomp, string double_ext) : Opencl_Module_Hmc() {
		Opencl_Module_Hmc::init(queue, 0, params, maxcomp, double_ext); /* init in body for proper this-pointer */
	};
	~Device() {
		finalize();
	};

  void runTestKernel(cl_mem in1, cl_mem in2, cl_mem out, cl_mem gf, int gs, int ls, int evenodd);
	void fill_kernels();
  void set_float_to_gm_squarenorm_device(cl_mem clmem_in, cl_mem clmem_out);
	void clear_kernels();
};

class Dummyfield : public Gaugefield_hybrid {

public:
	Dummyfield(cl_device_type device_type) : Gaugefield_hybrid() {
		std::stringstream tmp;
#ifdef _USEDOUBLEPREC_
		tmp << SOURCEDIR << "/tests/f_fermion_eo_input_1";
#else
		tmp << SOURCEDIR << "/tests/f_fermion_eo_input_1_single";
#endif
		params.readfile(tmp.str().c_str());

		init(1, device_type, &params);
	};

	virtual void init_tasks();
	virtual void finalize_opencl();

	hmc_float get_squarenorm(int which);
	void verify(hmc_float, hmc_float);	
	void runTestKernel(int evenodd);
        void reset_outfield();
private:
	void fill_buffers();
	void clear_buffers();
	inputparameters params;
  cl_mem in1, in2, in3, in4, out;
	cl_mem sqnorm;
	spinor * sf_in1;
  spinor * sf_in2;
	spinor * sf_in3;
  spinor * sf_in4;
	hmc_float * sf_out;
	Matrixsu3 * gf_in;

};

BOOST_AUTO_TEST_CASE( F_FERMION ){
	logger.info() << "Init CPU device";
	//params.print_info_inverter("m_gpu");
	Dummyfield cpu(CL_DEVICE_TYPE_CPU);
	logger.info() << "gaugeobservables: ";
	cpu.print_gaugeobservables_from_task(0, 0);
	logger.info() << "|phi_even_1|^2:";
	hmc_float cpu_back = cpu.get_squarenorm(0);
	logger.info() << "|phi_even_2|^2:";
	hmc_float cpu_back2 = cpu.get_squarenorm(1);
	logger.info() << "|phi_odd_1|^2:";
	hmc_float cpu_back3 = cpu.get_squarenorm(2);
	logger.info() << "|phi_odd_2|^2:";
	hmc_float cpu_back4 = cpu.get_squarenorm(3);
	cpu.runTestKernel(EVEN);
	logger.info() << "|force (even)|^2:";
	hmc_float cpu_res;
	cpu_res = cpu.get_squarenorm(4);
	cpu.reset_outfield();
	cpu.runTestKernel(ODD);
	logger.info() << "|force (odd)|^2:";
	hmc_float cpu_res2;
	cpu_res2 = cpu.get_squarenorm(4);
	cpu.runTestKernel(EVEN);
	logger.info() << "|force (even) + force (odd)|^2  (without setting the outvector to zero before):";
	hmc_float cpu_res3;
	cpu_res3 = cpu.get_squarenorm(4);
	BOOST_MESSAGE("Tested CPU");

	logger.info() << "Init GPU device";
	//params.print_info_inverter("m_gpu");
	Dummyfield gpu(CL_DEVICE_TYPE_GPU);
	logger.info() << "gaugeobservables: ";
	gpu.print_gaugeobservables_from_task(0, 0);
	logger.info() << "|phi_even_1|^2:";
	hmc_float gpu_back = gpu.get_squarenorm(0);
	logger.info() << "|phi_even_2|^2:";
	hmc_float gpu_back2 = gpu.get_squarenorm(1);
	logger.info() << "|phi_odd_1|^2:";
	hmc_float gpu_back3 = gpu.get_squarenorm(2);
	logger.info() << "|phi_odd_2|^2:";
	hmc_float gpu_back4 = gpu.get_squarenorm(3);
	gpu.runTestKernel(EVEN);
	logger.info() << "|force (even)|^2:";
	hmc_float gpu_res;
	gpu_res = gpu.get_squarenorm(4);
	gpu.reset_outfield();
	gpu.runTestKernel(ODD);
	logger.info() << "|force (odd)|^2:";
	hmc_float gpu_res2;
	gpu_res2 = gpu.get_squarenorm(4);
	gpu.runTestKernel(EVEN);
	logger.info() << "|force (even) + force (odd)|^2 (without setting the outvector to zero before):";
	hmc_float gpu_res3;
	gpu_res3 = gpu.get_squarenorm(4);
	BOOST_MESSAGE("Tested GPU");

	logger.info() << "Compare CPU and GPU results";
	logger.info() << "Input vectors:";
	cpu.verify(cpu_back, gpu_back);
	cpu.verify(cpu_back2, gpu_back2);
	cpu.verify(cpu_back3, gpu_back3);
	cpu.verify(cpu_back4, gpu_back4);
	logger.info() << "Output vectors:";
	cpu.verify(cpu_res, gpu_res);
	cpu.verify(cpu_res2, gpu_res2);
	cpu.verify(cpu_res3, gpu_res3);

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

void fill_sf_with_one(spinor * sf_in1, spinor * sf_in2, int size)
{
	for(int i = 0; i < size; ++i) {
		sf_in1[i].e0.e0 = hmc_complex_one;
		sf_in1[i].e0.e1 = hmc_complex_one;
		sf_in1[i].e0.e2 = hmc_complex_one;
		sf_in1[i].e1.e0 = hmc_complex_one;
		sf_in1[i].e1.e1 = hmc_complex_one;
		sf_in1[i].e1.e2 = hmc_complex_one;
		sf_in1[i].e2.e0 = hmc_complex_one;
		sf_in1[i].e2.e1 = hmc_complex_one;
		sf_in1[i].e2.e2 = hmc_complex_one;
		sf_in1[i].e3.e0 = hmc_complex_one;
		sf_in1[i].e3.e1 = hmc_complex_one;
		sf_in1[i].e3.e2 = hmc_complex_one;
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
{hmc_float val;
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

void fill_with_zero(hmc_float * sf_in, int size)
{
	for(int i = 0; i < size; ++i) {
	  sf_in[i] = 0.;
	}
	return;
}


void fill_sf_with_random(spinor * sf_in1, spinor * sf_in2, int size, int seed)
{
  //  Random rnd_loc(seed);
  for(int i = 0; i < size; ++i) {
  Random rnd_loc(seed);
			sf_in1[i].e0.e0.re = rnd_loc.doub();
			sf_in1[i].e0.e1.re = rnd_loc.doub();
			sf_in1[i].e0.e2.re = rnd_loc.doub();
			sf_in1[i].e1.e0.re = rnd_loc.doub();
			sf_in1[i].e1.e1.re = rnd_loc.doub();
			sf_in1[i].e1.e2.re = rnd_loc.doub();
			sf_in1[i].e2.e0.re = rnd_loc.doub();
			sf_in1[i].e2.e1.re = rnd_loc.doub();
			sf_in1[i].e2.e2.re = rnd_loc.doub();
			sf_in1[i].e3.e0.re = rnd_loc.doub();
			sf_in1[i].e3.e1.re = rnd_loc.doub();
			sf_in1[i].e3.e2.re = rnd_loc.doub();

			sf_in1[i].e0.e0.im = rnd_loc.doub();
			sf_in1[i].e0.e1.im = rnd_loc.doub();
			sf_in1[i].e0.e2.im = rnd_loc.doub();
			sf_in1[i].e1.e0.im = rnd_loc.doub();
			sf_in1[i].e1.e1.im = rnd_loc.doub();
			sf_in1[i].e1.e2.im = rnd_loc.doub();
			sf_in1[i].e2.e0.im = rnd_loc.doub();
			sf_in1[i].e2.e1.im = rnd_loc.doub();
			sf_in1[i].e2.e2.im = rnd_loc.doub();
			sf_in1[i].e3.e0.im = rnd_loc.doub();
			sf_in1[i].e3.e1.im = rnd_loc.doub();
			sf_in1[i].e3.e2.im = rnd_loc.doub();
			
			sf_in2[i] =sf_in1[i];

/*		  
			sf_in2[i].e0.e0.re = rnd_loc.doub();
			sf_in2[i].e0.e1.re = rnd_loc.doub();
			sf_in2[i].e0.e2.re = rnd_loc.doub();
			sf_in2[i].e1.e0.re = rnd_loc.doub();
			sf_in2[i].e1.e1.re = rnd_loc.doub();
			sf_in2[i].e1.e2.re = rnd_loc.doub();
			sf_in2[i].e2.e0.re = rnd_loc.doub();
			sf_in2[i].e2.e1.re = rnd_loc.doub();
			sf_in2[i].e2.e2.re = rnd_loc.doub();
			sf_in2[i].e3.e0.re = rnd_loc.doub();
			sf_in2[i].e3.e1.re = rnd_loc.doub();
			sf_in2[i].e3.e2.re = rnd_loc.doub();

			sf_in2[i].e0.e0.im = rnd_loc.doub();
			sf_in2[i].e0.e1.im = rnd_loc.doub();
			sf_in2[i].e0.e2.im = rnd_loc.doub();
			sf_in2[i].e1.e0.im = rnd_loc.doub();
			sf_in2[i].e1.e1.im = rnd_loc.doub();
			sf_in2[i].e1.e2.im = rnd_loc.doub();
			sf_in2[i].e2.e0.im = rnd_loc.doub();
			sf_in2[i].e2.e1.im = rnd_loc.doub();
			sf_in2[i].e2.e2.im = rnd_loc.doub();
			sf_in2[i].e3.e0.im = rnd_loc.doub();
			sf_in2[i].e3.e1.im = rnd_loc.doub();
			sf_in2[i].e3.e2.im = rnd_loc.doub();	       
*/
  }
  return;
}

void Dummyfield::reset_outfield(){
  size_t ae_buf_size = get_parameters()->get_gm_buf_size();
  int err = clEnqueueWriteBuffer(static_cast<Device*>(opencl_modules[0])->get_queue(), out, CL_TRUE, 0, ae_buf_size, sf_out, 0, 0, NULL);
  BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

  return;
}

void Dummyfield::fill_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	cl_int err;

	cl_context context = opencl_modules[0]->get_context();

	int NUM_ELEMENTS_SF;
	if(get_parameters()->get_use_eo() == true) NUM_ELEMENTS_SF =  params.get_eoprec_spinorfieldsize();
	else NUM_ELEMENTS_SF =  params.get_spinorfieldsize();

	int NUM_ELEMENTS_AE = params.get_gaugemomentasize()*params.get_su3algebrasize();

	sf_in1 = new spinor[NUM_ELEMENTS_SF];
	sf_in2 = new spinor[NUM_ELEMENTS_SF];
	sf_in3 = new spinor[NUM_ELEMENTS_SF];
	sf_in4 = new spinor[NUM_ELEMENTS_SF];
	sf_out = new hmc_float[NUM_ELEMENTS_AE];

	//use the variable use_cg to switch between cold and random input sf
	if(get_parameters()->get_use_cg() == true) {
	  fill_sf_with_one(sf_in1, sf_in2, NUM_ELEMENTS_SF);
	  fill_sf_with_one(sf_in3, sf_in4, NUM_ELEMENTS_SF);
	}
	else {
	  fill_sf_with_random(sf_in1, sf_in2, NUM_ELEMENTS_SF, 123456);
	  fill_sf_with_random(sf_in3, sf_in4, NUM_ELEMENTS_SF, 789101);
	}
	
	//filll with zeros
	hmc_float val = 0.;
	hmc_float val2 = 0.;
	fill_sf_with_float(sf_in1, NUM_ELEMENTS_SF, val);
	fill_sf_with_float(sf_in2, NUM_ELEMENTS_SF, val);
	fill_sf_with_float(sf_in3, NUM_ELEMENTS_SF, val2);
	fill_sf_with_float(sf_in4, NUM_ELEMENTS_SF, val2);
	

	fill_sf_with_pos(sf_in1,2/* NUM_ELEMENTS_SF*/);
	//fill_sf_with_pos(sf_in2, NUM_ELEMENTS_SF);
	fill_sf_with_pos(sf_in3, 2/*NUM_ELEMENTS_SF*/);
	//fill_sf_with_pos(sf_in4, NUM_ELEMENTS_SF);

	BOOST_REQUIRE(sf_in1);
	BOOST_REQUIRE(sf_in2);
	BOOST_REQUIRE(sf_in3);
	BOOST_REQUIRE(sf_in4);
	
	fill_with_zero(sf_out, NUM_ELEMENTS_AE);

	size_t sf_buf_size;
	if(get_parameters()->get_use_eo() == true) sf_buf_size = get_parameters()->get_eo_sf_buf_size();
	else sf_buf_size = get_parameters()->get_eo_sf_buf_size();
	size_t ae_buf_size = get_parameters()->get_gm_buf_size();
	//create buffer for sf on device (and copy sf_in to both for convenience)

	in1 = clCreateBuffer(context, CL_MEM_READ_ONLY , sf_buf_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	in2 = clCreateBuffer(context, CL_MEM_READ_ONLY , sf_buf_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	in3 = clCreateBuffer(context, CL_MEM_READ_ONLY , sf_buf_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	in4 = clCreateBuffer(context, CL_MEM_READ_ONLY , sf_buf_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	err = clEnqueueWriteBuffer(static_cast<Device*>(opencl_modules[0])->get_queue(), in1, CL_TRUE, 0, sf_buf_size, sf_in1, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	err = clEnqueueWriteBuffer(static_cast<Device*>(opencl_modules[0])->get_queue(), in2, CL_TRUE, 0, sf_buf_size, sf_in2, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	err = clEnqueueWriteBuffer(static_cast<Device*>(opencl_modules[0])->get_queue(), in3, CL_TRUE, 0, sf_buf_size, sf_in3, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	err = clEnqueueWriteBuffer(static_cast<Device*>(opencl_modules[0])->get_queue(), in4, CL_TRUE, 0, sf_buf_size, sf_in4, 0, 0, NULL);
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);

	out = clCreateBuffer(context, CL_MEM_WRITE_ONLY, ae_buf_size, 0, &err );
	BOOST_REQUIRE_EQUAL(err, CL_SUCCESS);
	this->reset_outfield();

	//create buffer for squarenorm on device
	sqnorm = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(hmc_float), 0, &err);
}

void Device::fill_kernels()
{
	//one only needs some kernels up to now. to save time during compiling they are put in here by hand
 	Opencl_Module::fill_kernels();

	//to this end, one has to set the needed files by hand
	basic_opencl_code = ClSourcePackage() << "opencl_header.cl" << "operations_geometry.cl" << "operations_complex.cl"
	                    << "operations_matrix_su3.cl" << "operations_matrix.cl" << "operations_gaugefield.cl";
	basic_fermion_code = basic_opencl_code << "types_fermions.h" << "operations_su3vec.cl" << "operations_spinor.cl" << "spinorfield.cl";
	//	basic_hmc_code = basic_fermion_code << "types_hmc.h";

	global_squarenorm_eoprec = createKernel("global_squarenorm_eoprec") << basic_fermion_code << "spinorfield_eo_squarenorm.cl";
	global_squarenorm_reduction = createKernel("global_squarenorm_reduction") << basic_fermion_code << "spinorfield_squarenorm.cl";

	ae_sqn = createKernel("gaugemomentum_squarenorm") << basic_fermion_code << "types_hmc.h" << "operations_gaugemomentum.cl" << "gaugemomentum_squarenorm.cl";

	testKernel = createKernel("fermion_force_eoprec") << basic_fermion_code << "types_hmc.h"  << "operations_gaugemomentum.cl" << "operations_spinorfield_eo.cl" << "fermionmatrix.cl" << "force_fermion_eo.cl";

}

void Dummyfield::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	clReleaseMemObject(in1);
	clReleaseMemObject(in2);
	clReleaseMemObject(in3);
	clReleaseMemObject(in4);
	clReleaseMemObject(out);
	clReleaseMemObject(sqnorm);

	delete[] sf_in1;
	delete[] sf_in2;
	delete[] sf_in3;
	delete[] sf_in4;
	delete[] sf_out;
}

void Device::clear_kernels()
{
	clReleaseKernel(testKernel);
	Opencl_Module::clear_kernels();
}

void Device::runTestKernel(cl_mem out, cl_mem in1, cl_mem in2, cl_mem gf, int gs, int ls, int evenodd)
{
	cl_int err;
	int eo;
	if(evenodd == ODD)
		eo=ODD;
	else
		eo=EVEN;
	err = clSetKernelArg(testKernel, 0, sizeof(cl_mem), &gf);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 1, sizeof(cl_mem), &in1);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 2, sizeof(cl_mem), &in2);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 3, sizeof(cl_mem), &out);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 4, sizeof(int), &eo);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	
	enqueueKernel(testKernel, gs, ls);
}

//this is a copy of "set_float_to_gaugemomentum_squarenorm_device"
void Device::set_float_to_gm_squarenorm_device(cl_mem clmem_in, cl_mem clmem_out)
{
  //__kernel void gaugemomentum_squarenorm(__global ae * in, __global hmc_float * out){
  //query work-sizes for kernel
  size_t ls2, gs2;
  cl_uint num_groups;
  this->get_work_sizes(ae_sqn, this->get_device_type(), &ls2, &gs2, &num_groups);
  //set arguments
  //__kernel void gaugemomentum_squarenorm(__global ae * in, __global hmc_float * out){
  int clerr = clSetKernelArg(ae_sqn, 0, sizeof(cl_mem), &clmem_in);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

  //  /** @todo add reduction */
  clerr = clSetKernelArg(ae_sqn,  1, sizeof(cl_mem), &clmem_out);
  if(clerr != CL_SUCCESS) throw Opencl_Error(clerr, "clSetKernelArg", __FILE__, __LINE__);

  enqueueKernel(ae_sqn  , gs2, ls2);
}


hmc_float Dummyfield::get_squarenorm(int which)
{
	//which controlls if the in or out-vector is looked at
	if(which == 0) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_eoprec_device(in1, sqnorm);
	if(which == 1) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_eoprec_device(in2, sqnorm);
	if(which == 2) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_eoprec_device(in3, sqnorm);
	if(which == 3) static_cast<Device*>(opencl_modules[0])->set_float_to_global_squarenorm_eoprec_device(in4, sqnorm);
	if(which == 4) static_cast<Device*>(opencl_modules[0])->set_float_to_gm_squarenorm_device(out, sqnorm);
	// get stuff from device
	hmc_float result;
	cl_int err = clEnqueueReadBuffer(*queue, sqnorm, CL_TRUE, 0, sizeof(hmc_float), &result, 0, 0, 0);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	logger.info() << result;
	return result;
}

void Dummyfield::verify(hmc_float cpu, hmc_float gpu){
  //this is too much required, since rounding errors can occur
  //  BOOST_REQUIRE_EQUAL(cpu, gpu);
  //instead, test if the two number agree within some percent
  hmc_float dev = (cpu - gpu)/cpu/100.;
  if(abs(dev) < 1e-10){
    logger.info() << "CPU and GPU result agree within accuary of " << 1e-10;
        logger.info() << "cpu: " << cpu << "\tgpu: " << gpu;
  }else{
    logger.info() << "CPU and GPU result DO NOT agree within accuary of " << 1e-10;
    logger.info() << "cpu: " << cpu << "\tgpu: " << gpu;
    BOOST_REQUIRE_EQUAL(1,0);
  }
}

void Dummyfield::runTestKernel(int evenodd)
{
	int gs, ls;
	if(opencl_modules[0]->get_device_type() == CL_DEVICE_TYPE_GPU) {
		gs = get_parameters()->get_eoprec_spinorfieldsize();
		ls = 64;
	} else if(opencl_modules[0]->get_device_type() == CL_DEVICE_TYPE_CPU) {
		gs = opencl_modules[0]->get_max_compute_units();
		ls = 1;
	}
	//interprete Y = (in1, in2) X = (in3, in4)
	//Y_odd = in2, Y_even = in1, X_odd = in4, X_even = in3
	if(evenodd==ODD){
		//this is then force(Y_odd, X_even) == force(in2, in3)
		static_cast<Device*>(opencl_modules[0])->runTestKernel(out, in2, in3, *(get_clmem_gaugefield()), gs, ls, evenodd);
	}
	else{
		//this is then force(Y_even, X_odd) == force(in1, in4)
		static_cast<Device*>(opencl_modules[0])->runTestKernel(out, in1, in4, *(get_clmem_gaugefield()), gs, ls, evenodd);
	}
}

