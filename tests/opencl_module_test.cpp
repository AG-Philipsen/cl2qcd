/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "kernelTester.hpp"

class PlaquetteTester : public KernelTester
{
public:
  PlaquetteTester(std::string inputfile, int typeOfPlaquette = 1):
    KernelTester("plaquette", inputfile), typeOfPlaquette(typeOfPlaquette)
  {
    callSpecificKernel();
  }
  void callSpecificKernel() override
  {
    auto device = this->system->get_devices()[0];
    auto * code = device->get_gaugefield_code(); 
    
    const hardware::buffers::Plain<hmc_float> plaq(1, device );
    const hardware::buffers::Plain<hmc_float> splaq(1, device);
    const hardware::buffers::Plain<hmc_float> tplaq(1, device);
    
    code->plaquette_device(this->gaugefield->get_buffers()[0], &plaq, &tplaq, &splaq);
	
    plaq.dump(&dev_plaq);
    splaq.dump(&dev_splaq);
    tplaq.dump(&dev_tplaq);
    
    switch( typeOfPlaquette )
      {
      case 1:
	kernelResult = dev_plaq;
	break;
      case 2:
	kernelResult = dev_tplaq;
	break;
      case 3:
	kernelResult = dev_splaq;
	break;
      default:
	throw std::invalid_argument(  "Do not recognize type of plaquette. Should be 1,2 or 3 (normal plaquette, temporal plaquette, spatial plaquette)" );
	break;
      }
  }
private:
  int typeOfPlaquette;
  hmc_float dev_plaq, dev_tplaq, dev_splaq;
};

class RectanglesTester : public KernelTester
{
public:
  RectanglesTester(std::string inputfile):
    KernelTester("rectangles", inputfile)
  {
    callSpecificKernel();
  }
  void callSpecificKernel() override
  {
    auto device = this->system->get_devices()[0];
    auto * code = device->get_gaugefield_code(); 
    
    hmc_float cpu_rect;
    code->gaugeobservables_rectangles(this->gaugefield->get_buffers()[0], &cpu_rect);
    kernelResult = cpu_rect;
  }
};

class StoutSmearTester : public KernelTester
{
public:
  StoutSmearTester(std::string inputfile):
    KernelTester("plaquette", inputfile)
  {
    callSpecificKernel();
  }
  void callSpecificKernel() override
  {
    auto device = this->system->get_devices()[0];
    auto * code = device->get_gaugefield_code(); 

    const hardware::buffers::Plain<hmc_float> plaq(1, device );
    const hardware::buffers::Plain<hmc_float> splaq(1, device);
    const hardware::buffers::Plain<hmc_float> tplaq(1, device);
    const hardware::buffers::SU3 out(this->gaugefield->get_buffers()[0]->get_elements(), device);

    code->stout_smear_device( this->gaugefield->get_buffers()[0], &out);

    code->plaquette_device( &out, &plaq, &tplaq, &splaq);
    plaq.dump(&dev_plaq);

    kernelResult = dev_plaq;
 }
private:
  hmc_float dev_plaq, dev_tplaq, dev_splaq;
  hmc_complex dev_pol;
};


#include "../meta/util.hpp"
#include "../physics/lattices/gaugefield.hpp"
#include "../hardware/device.hpp"

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE
#include <boost/test/unit_test.hpp>

#include "test_util.h"
class TestGaugefield {
public:
	TestGaugefield(const hardware::System * system) : system(system), prng(*system), gf(*system, prng) {
		BOOST_REQUIRE_EQUAL(system->get_devices().size(), 1);
		auto inputfile = system->get_inputparameters();
		std::string name = "test program";
		meta::print_info_hmc(inputfile);
		logger.info() << "gaugeobservables: ";
		print_gaugeobservables();
	};

	void print_gaugeobservables() {
		physics::lattices::print_gaugeobservables(gf, 0);
	}

	//todo: rename, this is just the buffer!
	const hardware::buffers::SU3 * get_gaugefield() {
		return gf.get_buffers()[0];
	}
	
	const meta::Inputparameters getParameters()
	{
		return system->get_inputparameters();
	}

	const hardware::code::Gaugefield * get_device();
	const hardware::System * const system;
	physics::PRNG prng;
	physics::lattices::Gaugefield gf;
};

const hardware::code::Gaugefield* TestGaugefield::get_device()
{
	return system->get_devices()[0]->get_gaugefield_code();
}

void gaugeobservables(const hardware::buffers::SU3 * gf, hmc_float * plaq_out, hmc_float * tplaq_out, hmc_float * splaq_out, hmc_complex * pol_out, const meta::Inputparameters& params)
{
	auto device = gf->get_device();
	auto code = device->get_gaugefield_code();

	const hardware::buffers::Plain<hmc_float> plaq(1, device);
	const hardware::buffers::Plain<hmc_float> splaq(1, device);
	const hardware::buffers::Plain<hmc_float> tplaq(1, device);
	const hardware::buffers::Plain<hmc_complex> pol(1, device);

	//measure plaquette
	code->plaquette_device(gf, &plaq, &tplaq, &splaq);

	//read out values
	hmc_float tmp_plaq = 0.;
	hmc_float tmp_splaq = 0.;
	hmc_float tmp_tplaq = 0.;
	//NOTE: these are blocking calls!
	plaq.dump(&tmp_plaq);
	splaq.dump(&tmp_splaq);
	tplaq.dump(&tmp_tplaq);

	tmp_tplaq /= static_cast<hmc_float>(meta::get_tplaq_norm(params));
	tmp_splaq /= static_cast<hmc_float>(meta::get_splaq_norm(params));
	tmp_plaq  /= static_cast<hmc_float>(meta::get_plaq_norm(params));

	(*plaq_out) = tmp_plaq;
	(*splaq_out) = tmp_splaq;
	(*tplaq_out) = tmp_tplaq;

	//measure polyakovloop
	code->polyakov_device(gf, &pol);

	//read out values
	hmc_complex tmp_pol = hmc_complex_zero;
	//NOTE: this is a blocking call!
	pol.dump(&tmp_pol);

	tmp_pol.re /= static_cast<hmc_float> ( meta::get_poly_norm(params) );
	tmp_pol.im /= static_cast<hmc_float> ( meta::get_poly_norm(params) );

	pol_out->re = tmp_pol.re;
	pol_out->im = tmp_pol.im;
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
	gaugeobservables(dummy.get_gaugefield(), &dev_plaq, &dev_tplaq, &dev_splaq, &dev_pol, params);

	logger.info() << "reference value:\t" << "values obtained from host functionality";
	hmc_float prec = params.get_solver_prec();
	logger.info() << "acceptance precision: " << prec;

	// verify
	BOOST_CHECK_CLOSE(ref_pol.re, dev_pol.re,   prec);
	BOOST_CHECK_CLOSE(ref_pol.im, dev_pol.im,   prec);
}

BOOST_AUTO_TEST_SUITE ( PLAQUETTE )

BOOST_AUTO_TEST_CASE( PLAQUETTE_1 )
{
  PlaquetteTester plaquetteTester("plaquette_input_1", 1);
}

BOOST_AUTO_TEST_CASE( PLAQUETTE_2 )
{
  PlaquetteTester plaquetteTester("plaquette_input_2", 1);
}

BOOST_AUTO_TEST_CASE( PLAQUETTE_3 )
{
  PlaquetteTester plaquetteTester("plaquette_input_3", 1);
}

BOOST_AUTO_TEST_CASE( PLAQUETTE_TEMPORAL_1 )
{
  PlaquetteTester plaquetteTester("plaquette_temporal_input_1", 2);
}

BOOST_AUTO_TEST_CASE( PLAQUETTE_TEMPORAL_2 )
{
  PlaquetteTester plaquetteTester("plaquette_temporal_input_2", 2);
}

BOOST_AUTO_TEST_CASE( PLAQUETTE_TEMPORAL_3 )
{
  PlaquetteTester plaquetteTester("plaquette_temporal_input_3", 2);
}

BOOST_AUTO_TEST_CASE( PLAQUETTE_SPATIAL_1 )
{
  PlaquetteTester plaquetteTester("plaquette_spatial_input_1", 3);
}

BOOST_AUTO_TEST_CASE( PLAQUETTE_SPATIAL_2 )
{
  PlaquetteTester plaquetteTester("plaquette_spatial_input_2", 3);
}

BOOST_AUTO_TEST_CASE( PLAQUETTE_SPATIAL_3 )
{
  PlaquetteTester plaquetteTester("plaquette_spatial_input_3", 3);
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( PLAQUETTE_REDUCTION )

BOOST_AUTO_TEST_CASE( PLAQUETTE_REDUCTION_1 )
{
  PlaquetteTester plaquetteTester("plaquette_reduction_input_1");
}

BOOST_AUTO_TEST_CASE( PLAQUETTE_REDUCTION_2 )
{
  PlaquetteTester plaquetteTester("plaquette_reduction_input_2");
}

BOOST_AUTO_TEST_CASE( PLAQUETTE_REDUCTION_3 )
{
  PlaquetteTester plaquetteTester("plaquette_reduction_input_3");
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

BOOST_AUTO_TEST_SUITE ( POLYAKOV_REDUCTION )

BOOST_AUTO_TEST_CASE( POLYAKOV_REDUCTION_1 )
{
	test_polyakov( "/polyakov_reduction_input_1", {1, 0} );
}

BOOST_AUTO_TEST_CASE( POLYAKOV_REDUCTION_2 )
{
	test_polyakov( "/polyakov_reduction_input_2", {1, 0} );
}

BOOST_AUTO_TEST_CASE( POLYAKOV_REDUCTION_3 )
{
	test_polyakov( "/polyakov_reduction_input_3", {1, 0} );
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
  StoutSmearTester StoutSmearTester("stout_smear_input_1");
}

BOOST_AUTO_TEST_CASE( STOUT_SMEAR_2 )
{
  StoutSmearTester StoutSmearTester("stout_smear_input_2");
}

BOOST_AUTO_TEST_CASE( STOUT_SMEAR_3 )
{
  StoutSmearTester StoutSmearTester("stout_smear_input_3");
}

BOOST_AUTO_TEST_CASE( STOUT_SMEAR_4 )
{
  StoutSmearTester StoutSmearTester("stout_smear_input_4");
}

BOOST_AUTO_TEST_CASE( STOUT_SMEAR_5 )
{
  StoutSmearTester StoutSmearTester("stout_smear_input_5");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( RECTANGLES )

BOOST_AUTO_TEST_CASE( RECTANGLES_1 )
{
  RectanglesTester rectanglesTester("rectangles_input_1");
}

BOOST_AUTO_TEST_CASE( RECTANGLES_2 )
{
  RectanglesTester rectanglesTester("rectangles_input_2");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( RECTANGLES_REDUCTION )

BOOST_AUTO_TEST_CASE( RECTANGLES_REDUCTION_1 )
{
  RectanglesTester rectanglesTester("rectangles_reduction_input_1");
}

BOOST_AUTO_TEST_CASE( RECTANGLES_REDUCTION_2 )
{
  RectanglesTester rectanglesTester("rectangles_reduction_input_2");
}

BOOST_AUTO_TEST_CASE( RECTANGLES_REDUCTION_3 )
{
  RectanglesTester rectanglesTester("rectangles_reduction_input_3");
}

BOOST_AUTO_TEST_SUITE_END()

