#include "../gaugefield_hybrid.h"
#include "../meta/util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE stout_smear_test
#include <boost/test/unit_test.hpp>

//some functionality                                                                                             
#include "test_util.h"

std::string const exec_name = "stout_smear_test";

#define CLX_CHECK_CLOSE(left, right, precision) \
{ \
  BOOST_CHECK_CLOSE(left.re, right.re, precision); \
  BOOST_CHECK_CLOSE(left.im, right.im, precision); \
}

//global parameters
hmc_float rho = 0.01;
int iter = 1;

class Device : public hardware::code::Opencl_Module {
protected:
	virtual size_t get_read_write_size(const std::string&) const {
		return 0;
	};
	virtual uint64_t get_flop_size(const std::string&) const {
		return 0;
	};
public:
	Device(const meta::Inputparameters& params, hardware::Device * device) : Opencl_Module(params, device), out(get_device()->get_gaugefield_code()->get_gaugefield()->get_elements(), device) {
	};

	const hardware::buffers::SU3 out;
};

class Dummyfield : public Gaugefield_hybrid {

public:
	Dummyfield(const hardware::System * system) : Gaugefield_hybrid(system) {
		auto inputfile = system->get_inputparameters();
		init(1, inputfile.get_use_gpu() ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU);
		meta::print_info_hmc(exec_name.c_str(), inputfile);


		logger.info() << "gaugeobservables: ";
		this->print_gaugeobservables_from_task(0, 0);



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
	opencl_modules = new hardware::code::Opencl_Module* [get_num_tasks()];
	opencl_modules[0] = new Device(get_parameters(), get_device_for_task(0));

	fill_buffers();
}

void Dummyfield::finalize_opencl()
{
	Gaugefield_hybrid::finalize_opencl();
}

void Dummyfield::fill_buffers()
{
  //set_gaugefield_cold(get_sgf());
	copy_gaugefield_to_task(0);
}

void Dummyfield::runTestKernel()
{
	Device* device = static_cast<Device*>(opencl_modules[0]);
	//CP: this currently causes a segfault!!!
	auto gf_code = get_device_for_task(0)->get_gaugefield_code();
	gf_code->stout_smear_device( gf_code->get_gaugefield(), &device->out);
}

void Dummyfield::get_gaugeobservables_from_task(int ntask, hmc_float * plaq, hmc_float * tplaq, hmc_float * splaq, hmc_complex * pol)
{
	if( ntask < 0 || ntask > get_num_tasks() ) throw Print_Error_Message("devicetypes index out of range", __FILE__, __LINE__);
	auto gf_code = get_device_for_task(ntask)->get_gaugefield_code();
	gf_code->gaugeobservables(plaq, tplaq, splaq, pol);
}

//this is just out of laziness, a copy of the function above
void Dummyfield::get_gaugeobservables_from_task(int dummy, int ntask, hmc_float * plaq, hmc_float * tplaq, hmc_float * splaq, hmc_complex * pol)
{
	dummy = 0;
	if( ntask < 0 || ntask > get_num_tasks() ) throw Print_Error_Message("devicetypes index out of range", __FILE__, __LINE__);
	auto gf_code = get_device_for_task(dummy)->get_gaugefield_code();
	gf_code->gaugeobservables(&static_cast<Device*>(opencl_modules[0])->out, plaq, tplaq, splaq, pol);
}

class TestGaugefield : public Gaugefield_hybrid {
public:
	TestGaugefield(const hardware::System * system) : Gaugefield_hybrid(system) {
		auto inputfile = system->get_inputparameters();
		init(1, inputfile.get_use_gpu() ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU);
		std::string name = "test program";
		meta::print_info_hmc(name.c_str(), inputfile);
		logger.info() << "gaugeobservables: ";
		this->print_gaugeobservables_from_task(0, 0);
	};
	virtual void finalize_opencl() override;

	hardware::code::Gaugefield * get_device();
};

hardware::code::Gaugefield* TestGaugefield::get_device()
{
	return static_cast<hardware::code::Gaugefield*>(opencl_modules[0]);
}

void TestGaugefield::finalize_opencl()
{
	Gaugefield_hybrid::finalize_opencl();
}

BOOST_AUTO_TEST_CASE( STOUT_SMEAR )
{
  std::string inputfile = "/stout_smear_input_1";
	std::string kernelName = "stout_smear";
        printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	//Dummyfield dummy(&system);

	TestGaugefield dummy2(&system);
	hardware::code::Gaugefield * device = dummy2.get_device();
	auto gf_code = device->get_device()->get_gaugefield_code();
	//out buffer
	const hardware::buffers::SU3 out(device->get_gaugefield()->get_elements(), device->get_device());
	hmc_float plaq_cpu, tplaq_cpu, splaq_cpu;
	hmc_complex pol_cpu;

	logger.info() << "gaugeobservables of in field before: ";
	dummy2.print_gaugeobservables_from_task(0,0);
	//dummy.print_gaugeobservables_from_task(0, 0);
	logger.info() << "gaugeobservables of out field before: ";
	//dummy.get_gaugeobservables_from_task(0, 0, &plaq_cpu, &tplaq_cpu, &splaq_cpu, &pol_cpu);
	//logger.info() << "plaq: " << plaq_cpu << "\t" << tplaq_cpu  << "\t" << splaq_cpu  << "\t" << pol_cpu.re  << "\t" << pol_cpu.im ;
        gf_code->gaugeobservables(&out , &plaq_cpu, &tplaq_cpu, &splaq_cpu, &pol_cpu);
	logger.info() << "plaq: " << plaq_cpu << "\t" << tplaq_cpu  << "\t" << splaq_cpu  << "\t" << pol_cpu.re  << "\t" << pol_cpu.im ;

	//dummy.runTestKernel();
	gf_code->stout_smear_device( gf_code->get_gaugefield(), &out);

	logger.info() << "gaugeobservables of in field after: ";
	//dummy.print_gaugeobservables_from_task(0, 0);
	dummy2.print_gaugeobservables_from_task(0, 0);
	logger.info() << "gaugeobservables of out field after: ";
	//dummy.get_gaugeobservables_from_task(0, 0, &plaq_cpu, &tplaq_cpu, &splaq_cpu, &pol_cpu);
	//logger.info() << "plaq: " << plaq_cpu << "\t" << tplaq_cpu  << "\t" << splaq_cpu  << "\t" << pol_cpu.re  << "\t" << pol_cpu.im ;

        gf_code->gaugeobservables(&out , &plaq_cpu, &tplaq_cpu, &splaq_cpu, &pol_cpu);
	logger.info() << "plaq: " << plaq_cpu << "\t" << tplaq_cpu  << "\t" << splaq_cpu  << "\t" << pol_cpu.re  << "\t" << pol_cpu.im ;

	/*
	logger.info() << "Choosing reference value and acceptance precision";
	hmc_float ref_val = params.get_test_ref_value();
	logger.info() << "reference value:\t" << ref_val;
	hmc_float prec = params.get_solver_prec();
	logger.info() << "acceptance precision: " << prec;
	
	logger.info() << "Compare result to reference value";
	BOOST_REQUIRE_CLOSE(plaq_cpu, ref_val, prec);
	logger.info() << "Done";
	*/
        testFloatAgainstInputparameters(plaq_cpu, params);
        BOOST_MESSAGE("Test done");
	BOOST_MESSAGE("Test done");
}
