#include "../gaugefield_hybrid.h"
#include "../meta/util.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE
#include <boost/test/unit_test.hpp>

//some functionality
#include "test_util.h"

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


void test_rectangles(std::string inputfile)
{
	std::string kernelName = "rectangles";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield cpu(&system);
	hardware::code::Gaugefield * device = cpu.get_device();

	logger.info() << "calc rectangles value:";
	hmc_float cpu_rect;
	device->gaugeobservables_rectangles(device->get_gaugefield(), &cpu_rect);
	logger.info() << cpu_rect;

	logger.info() << "Finalize device";
	cpu.finalize();

	testFloatAgainstInputparameters(cpu_rect, params);
	BOOST_MESSAGE("Test done");
}

void test_plaquette(std::string inputfile)
{
	std::string kernelName = "plaquette";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield dummy(&system);

	// get reference solutions
	hmc_float ref_plaq, ref_tplaq, ref_splaq;
	ref_plaq = dummy.plaquette(&ref_tplaq, &ref_splaq);

	// get device colutions
	hmc_float dev_plaq, dev_tplaq, dev_splaq;
	hmc_complex dev_pol;
	hardware::code::Gaugefield * device = dummy.get_device();
	device->gaugeobservables(&dev_plaq, &dev_tplaq, &dev_splaq, &dev_pol);

	logger.info() << "Finalize device";
	dummy.finalize();

	logger.info() << "reference value:\t" << "values obtained from host functionality";
	hmc_float prec = params.get_solver_prec();
	logger.info() << "acceptance precision: " << prec;

	// verify
	BOOST_REQUIRE_CLOSE(ref_plaq,   dev_plaq,     prec);
	BOOST_REQUIRE_CLOSE(ref_tplaq,  dev_tplaq,    prec);
	BOOST_REQUIRE_CLOSE(ref_splaq,  dev_splaq,    prec);
}

void test_polyakov(std::string inputfile)
{
	std::string kernelName = "plaquette";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield dummy(&system);

	// get reference solutions
	hmc_complex ref_pol = dummy.polyakov();

	// get device colutions
	hmc_float dev_plaq, dev_tplaq, dev_splaq;
	hmc_complex dev_pol;
	hardware::code::Gaugefield * device = dummy.get_device();
	device->gaugeobservables(&dev_plaq, &dev_tplaq, &dev_splaq, &dev_pol);

	logger.info() << "Finalize device";
	dummy.finalize();

	logger.info() << "reference value:\t" << "values obtained from host functionality";
	hmc_float prec = params.get_solver_prec();
	logger.info() << "acceptance precision: " << prec;

	// verify
	BOOST_REQUIRE_CLOSE(ref_pol.re, dev_pol.re,   prec);
	BOOST_REQUIRE_CLOSE(ref_pol.im, dev_pol.im,   prec);
}

void test_stout_smear(std::string inputfile)
{
  std::string kernelName = "stout_smear";
  printKernelInfo(kernelName);
  logger.info() << "Init device";
  meta::Inputparameters params = create_parameters(inputfile);
  hardware::System system(params);

  TestGaugefield dummy2(&system);
  hardware::code::Gaugefield * device = dummy2.get_device();
  auto gf_code = device->get_device()->get_gaugefield_code();
  //out buffer                                                                                                                                                                
  const hardware::buffers::SU3 out(device->get_gaugefield()->get_elements(), device->get_device());
  hmc_float plaq_cpu, tplaq_cpu, splaq_cpu;
  hmc_complex pol_cpu;

  logger.info() << "gaugeobservables of in field before: ";
  dummy2.print_gaugeobservables_from_task(0,0);
  logger.info() << "gaugeobservables of out field before: ";
  gf_code->gaugeobservables(&out , &plaq_cpu, &tplaq_cpu, &splaq_cpu, &pol_cpu);
  logger.info() << "plaq: " << plaq_cpu << "\t" << tplaq_cpu  << "\t" << splaq_cpu  << "\t" << pol_cpu.re  << "\t" << pol_cpu.im ;

  gf_code->stout_smear_device( gf_code->get_gaugefield(), &out);

  logger.info() << "gaugeobservables of in field after: ";
  dummy2.print_gaugeobservables_from_task(0, 0);
  logger.info() << "gaugeobservables of out field after: ";

  gf_code->gaugeobservables(&out , &plaq_cpu, &tplaq_cpu, &splaq_cpu, &pol_cpu);
  logger.info() << "plaq: " << plaq_cpu << "\t" << tplaq_cpu  << "\t" << splaq_cpu  << "\t" << pol_cpu.re  << "\t" << pol_cpu.im ;

  testFloatAgainstInputparameters(plaq_cpu, params);
  BOOST_MESSAGE("Test done");

}

BOOST_AUTO_TEST_SUITE ( PLAQUETTE )

BOOST_AUTO_TEST_CASE( PLAQUETTE_1 )
{
	test_plaquette( "/plaquette_input_1" );
}

BOOST_AUTO_TEST_CASE( PLAQUETTE_2 )
{
	test_plaquette( "/plaquette_input_2" );
}

BOOST_AUTO_TEST_CASE( PLAQUETTE_3 )
{
	test_plaquette( "/plaquette_input_3" );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( POLYAKOV )

BOOST_AUTO_TEST_CASE( POLYAKOV_1 )
{
	test_polyakov( "/polyakov_input_1" );
}

BOOST_AUTO_TEST_CASE( POLYAKOV_2)
{
	test_polyakov( "/polyakov_input_2");
}

BOOST_AUTO_TEST_CASE( POLYAKOV_3)
{
	test_polyakov( "/polyakov_input_3");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( CONVERT_GAUGEFIELD_TO_SOA )

BOOST_AUTO_TEST_CASE( CONVERT_GAUGEFIELD_TO_SOA_1 )
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( CONVERT_GAUGEFIELD_FROM_SOA )

BOOST_AUTO_TEST_CASE( CONVERT_GAUGEFIELD_FROM_SOA_1 )
{
	BOOST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( STOUT_SMEAR )

BOOST_AUTO_TEST_CASE( STOUT_SMEAR_1 )
{
        test_stout_smear("/stout_smear_input_1");
}

BOOST_AUTO_TEST_CASE( STOUT_SMEAR_2 )
{
        test_stout_smear("/stout_smear_input_2");
}

BOOST_AUTO_TEST_CASE( STOUT_SMEAR_3 )
{
        test_stout_smear("/stout_smear_input_3");
}

BOOST_AUTO_TEST_CASE( STOUT_SMEAR_4 )
{
        test_stout_smear("/stout_smear_input_4");
}

BOOST_AUTO_TEST_CASE( STOUT_SMEAR_5 )
{
        test_stout_smear("/stout_smear_input_5");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( RECTANGLES )

BOOST_AUTO_TEST_CASE( RECTANGLES_1 )
{
	test_rectangles("/rectangles_input_1");
}

BOOST_AUTO_TEST_CASE( RECTANGLES_2 )
{
	test_rectangles("/rectangles_input_2");
}

BOOST_AUTO_TEST_SUITE_END()


