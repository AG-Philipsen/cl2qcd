#include "../meta/util.hpp"
#include "../physics/lattices/gaugefield.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE
#include <boost/test/unit_test.hpp>

//some functionality
#include "test_util.h"

class TestGaugefield {
public:
	TestGaugefield(const hardware::System * system) : system(system), prng(*system), gf(*system, prng) {
		BOOST_REQUIRE_EQUAL(system->get_devices().size(), 1);
		auto inputfile = system->get_inputparameters();
		std::string name = "test program";
		meta::print_info_hmc(name.c_str(), inputfile);
		logger.info() << "gaugeobservables: ";
		print_gaugeobservables();
	};

	void print_gaugeobservables() {
		physics::lattices::print_gaugeobservables(gf, 0);
	}

	const hardware::buffers::SU3 * get_gaugefield() {
		return gf.get_buffers()[0];
	}

	hardware::code::Gaugefield * get_device();
private:
	const hardware::System * const system;
	physics::PRNG prng;
	physics::lattices::Gaugefield gf;
};

hardware::code::Gaugefield* TestGaugefield::get_device()
{
	return system->get_devices()[0]->get_gaugefield_code();
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
	device->gaugeobservables_rectangles(cpu.get_gaugefield(), &cpu_rect);
	logger.info() << cpu_rect;

	testFloatAgainstInputparameters(cpu_rect, params);
	BOOST_MESSAGE("Test done");
}

void test_plaquette(std::string inputfile, hmc_float ref_plaq, hmc_float ref_tplaq, hmc_float ref_splaq)
{
	std::string kernelName = "plaquette";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield dummy(&system);

	// get device colutions
	hmc_float dev_plaq, dev_tplaq, dev_splaq;
	hmc_complex dev_pol;
	hardware::code::Gaugefield * device = dummy.get_device();
	device->gaugeobservables(dummy.get_gaugefield(), &dev_plaq, &dev_tplaq, &dev_splaq, &dev_pol);

	logger.info() << "reference value:\t" << "values obtained from host functionality";
	hmc_float prec = params.get_solver_prec();
	logger.info() << "acceptance precision: " << prec;

	// verify
	BOOST_CHECK_CLOSE(ref_plaq,   dev_plaq,     prec);
	BOOST_CHECK_CLOSE(ref_tplaq,  dev_tplaq,    prec);
	BOOST_CHECK_CLOSE(ref_splaq,  dev_splaq,    prec);
}

void test_polyakov(std::string inputfile, hmc_complex ref_pol)
{
	std::string kernelName = "plaquette";
	printKernelInfo(kernelName);
	logger.info() << "Init device";
	meta::Inputparameters params = create_parameters(inputfile);
	hardware::System system(params);
	TestGaugefield dummy(&system);

	// get device colutions
	hmc_float dev_plaq, dev_tplaq, dev_splaq;
	hmc_complex dev_pol;
	hardware::code::Gaugefield * device = dummy.get_device();
	device->gaugeobservables(dummy.get_gaugefield(), &dev_plaq, &dev_tplaq, &dev_splaq, &dev_pol);

	logger.info() << "reference value:\t" << "values obtained from host functionality";
	hmc_float prec = params.get_solver_prec();
	logger.info() << "acceptance precision: " << prec;

	// verify
	BOOST_CHECK_CLOSE(ref_pol.re, dev_pol.re,   prec);
	BOOST_CHECK_CLOSE(ref_pol.im, dev_pol.im,   prec);
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
	const hardware::buffers::SU3 out(dummy2.get_gaugefield()->get_elements(), device->get_device());
	hmc_float plaq_cpu, tplaq_cpu, splaq_cpu;
	hmc_complex pol_cpu;

	logger.info() << "gaugeobservables of in field before: ";
	dummy2.print_gaugeobservables();
	logger.info() << "gaugeobservables of out field before: ";
	gf_code->gaugeobservables(&out , &plaq_cpu, &tplaq_cpu, &splaq_cpu, &pol_cpu);
	logger.info() << "plaq: " << plaq_cpu << "\t" << tplaq_cpu  << "\t" << splaq_cpu  << "\t" << pol_cpu.re  << "\t" << pol_cpu.im ;

	gf_code->stout_smear_device( dummy2.get_gaugefield(), &out);

	logger.info() << "gaugeobservables of in field after: ";
	dummy2.print_gaugeobservables();
	logger.info() << "gaugeobservables of out field after: ";

	gf_code->gaugeobservables(&out , &plaq_cpu, &tplaq_cpu, &splaq_cpu, &pol_cpu);
	logger.info() << "plaq: " << plaq_cpu << "\t" << tplaq_cpu  << "\t" << splaq_cpu  << "\t" << pol_cpu.re  << "\t" << pol_cpu.im ;

	testFloatAgainstInputparameters(plaq_cpu, params);
	BOOST_MESSAGE("Test done");

}

BOOST_AUTO_TEST_SUITE ( PLAQUETTE )

BOOST_AUTO_TEST_CASE( PLAQUETTE_1 )
{
	test_plaquette( "/plaquette_input_1", 1, 1, 1 );
}

BOOST_AUTO_TEST_CASE( PLAQUETTE_2 )
{
	test_plaquette( "/plaquette_input_2", 0.0050057845805392071, 0.00096087997130853749, 0.0090506891897698793 );
}

BOOST_AUTO_TEST_CASE( PLAQUETTE_3 )
{
	test_plaquette( "/plaquette_input_3", 0.57107711169452691, 0.57147433845588391, 0.57067988493316968 );
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( POLYAKOV )

BOOST_AUTO_TEST_CASE( POLYAKOV_1 )
{
	test_polyakov( "/polyakov_input_1", {1, 0} );
}

BOOST_AUTO_TEST_CASE( POLYAKOV_2)
{
	test_polyakov( "/polyakov_input_2", {0.01190176915346167, -0.073648068466866529} );
}

BOOST_AUTO_TEST_CASE( POLYAKOV_3)
{
	test_polyakov( "/polyakov_input_3", { -0.11349672123636854, 0.22828243566855222} );
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


