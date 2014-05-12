/*
 * Copyright 2012, 2013, 2014 Christopher Pinke, Matthias Bach
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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE hardware::code::Gaugefield
#include <boost/test/unit_test.hpp>

#include "kernelTester.hpp"

class GaugefieldTesterDouble : public KernelTesterDouble
{
public:
  GaugefieldTesterDouble(std::string kernelName, std::string inputfile):
    KernelTesterDouble(kernelName, inputfile)
  {
    prng = new physics::PRNG(*system);
    gaugefield = new physics::lattices::Gaugefield(*system, *prng);

    code = device->get_gaugefield_code();
  }
protected:
  const hardware::buffers::SU3* getGaugefieldBuffer()
  {
    return gaugefield->get_buffers()[0];
  }

  const hardware::code::Gaugefield * code;  

  physics::PRNG * prng;  
  physics::lattices::Gaugefield * gaugefield;
};

class PlaquetteTester : public GaugefieldTesterDouble
{
public:
  PlaquetteTester(std::string inputfile, int typeOfPlaquette = 1):
    GaugefieldTesterDouble("plaquette", inputfile), typeOfPlaquette(typeOfPlaquette)
  {
    const hardware::buffers::Plain<hmc_float> plaq(1, device );
    const hardware::buffers::Plain<hmc_float> splaq(1, device);
    const hardware::buffers::Plain<hmc_float> tplaq(1, device);
    
    code->plaquette_device(getGaugefieldBuffer(), &plaq, &tplaq, &splaq);
	
    switch( typeOfPlaquette )
      {
      case 1:
	plaq.dump(&kernelResult);
	break;
      case 2:
	tplaq.dump(&kernelResult);
	break;
      case 3:
	splaq.dump(&kernelResult);
	break;
      default:
	throw std::invalid_argument(  "Do not recognize type of plaquette. Should be 1,2 or 3 (normal plaquette, temporal plaquette, spatial plaquette)" );
	break;
      }
  }
private:
  int typeOfPlaquette;
};

class RectanglesTester : public GaugefieldTesterDouble
{
public:
  RectanglesTester(std::string inputfile):
    GaugefieldTesterDouble("rectangles", inputfile)
  {
    const hardware::buffers::Plain<hmc_float> rect(1, device );
    code->rectangles_device(getGaugefieldBuffer(), &rect);
    rect.dump(&kernelResult);
  }
};

class StoutSmearTester : public GaugefieldTesterDouble
{
public:
  StoutSmearTester(std::string inputfile):
    GaugefieldTesterDouble("stout_smear", inputfile)
  {
    auto gaugefieldBuffer = getGaugefieldBuffer();

    const hardware::buffers::Plain<hmc_float> plaq(1, device );
    const hardware::buffers::Plain<hmc_float> splaq(1, device);
    const hardware::buffers::Plain<hmc_float> tplaq(1, device);
    const hardware::buffers::SU3 out(gaugefieldBuffer->get_elements(), device);

    code->stout_smear_device( gaugefieldBuffer, &out);

    code->plaquette_device( &out, &plaq, &tplaq, &splaq);
    plaq.dump(&kernelResult);
 }
};

class PolyakovloopTester : public KernelTesterComplex
{
public:
  PolyakovloopTester(std::string inputfile):
    KernelTesterComplex("polyakov", inputfile)
  {
    prng = new physics::PRNG(*system);
    gaugefield = new physics::lattices::Gaugefield(*system, *prng);

    auto * code = device->get_gaugefield_code(); 

    const hardware::buffers::Plain<hmc_complex> pol(1, device);
    code->polyakov_device(getGaugefieldBuffer(), &pol);
    pol.dump(&kernelResult);
  }
private:
  physics::PRNG * prng;  
  physics::lattices::Gaugefield * gaugefield;

  const hardware::buffers::SU3* getGaugefieldBuffer()
  {
    return gaugefield->get_buffers()[0];
  }
};

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
  PolyakovloopTester polyakovloopTester("/polyakov_input_1");
}

BOOST_AUTO_TEST_CASE( POLYAKOV_2)
{
  PolyakovloopTester polyakovloopTester("/polyakov_input_2");
}

BOOST_AUTO_TEST_CASE( POLYAKOV_3)
{
  PolyakovloopTester polyakovloopTester("/polyakov_input_3");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( POLYAKOV_REDUCTION )

BOOST_AUTO_TEST_CASE( POLYAKOV_REDUCTION_1 )
{
  PolyakovloopTester polyakovloopTester("/polyakov_reduction_input_1");
}

BOOST_AUTO_TEST_CASE( POLYAKOV_REDUCTION_2 )
{
  PolyakovloopTester polyakovloopTester("/polyakov_reduction_input_2");
}

BOOST_AUTO_TEST_CASE( POLYAKOV_REDUCTION_3 )
{
  PolyakovloopTester polyakovloopTester("/polyakov_reduction_input_3");
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

