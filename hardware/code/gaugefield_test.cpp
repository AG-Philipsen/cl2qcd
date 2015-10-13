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
#include "../../physics/prng.hpp"
#include "../../physics/lattices/gaugefield.hpp"
#include "../buffers/plain.hpp"
#include "gaugefield.hpp"

class GaugefieldTester : public KernelTester {
public:
	GaugefieldTester(std::string kernelName, std::string inputfileIn, int numberOfValues = 1):
	  KernelTester(kernelName, getSpecificInputfile(inputfileIn), numberOfValues), prngParameters(&system->get_inputparameters()) {
		prng = new physics::PRNG(*system, &prngParameters);
		params = new physics::lattices::GaugefieldParametersImplementation( &system->get_inputparameters());
		gaugefield = new physics::lattices::Gaugefield(*system, params, *prng);
		code = device->getGaugefieldCode();
	}

	~GaugefieldTester()
	{
		delete prng;
		delete params;
		delete gaugefield;
	}

protected:
	std::string getSpecificInputfile(std::string inputfileIn)
	{
		return "gaugefield/" + inputfileIn;
	}

	const hardware::buffers::SU3* getGaugefieldBuffer() {
		return gaugefield->get_buffers()[0];
	}

	physics::PRNG * prng;
	const hardware::code::Gaugefield * code;
	physics::lattices::Gaugefield * gaugefield;
	physics::lattices::GaugefieldParametersImplementation * params;
	physics::ParametersPrng_fromMetaInputparameters prngParameters;
};

BOOST_AUTO_TEST_SUITE ( PLAQUETTE )

	class PlaquetteTester : public GaugefieldTester {
	public:
		PlaquetteTester(std::string inputfile, int typeOfPlaquette = 1):
			GaugefieldTester("plaquette", inputfile, 1), typeOfPlaquette(typeOfPlaquette) {
			const hardware::buffers::Plain<hmc_float> plaq(1, device );
			const hardware::buffers::Plain<hmc_float> splaq(1, device);
			const hardware::buffers::Plain<hmc_float> tplaq(1, device);

			code->plaquette_device(getGaugefieldBuffer(), &plaq, &tplaq, &splaq);

			switch( typeOfPlaquette ) {
				case 1:
					plaq.dump(&kernelResult[0]);
					break;
				case 2:
					tplaq.dump(&kernelResult[0]);
					break;
				case 3:
					splaq.dump(&kernelResult[0]);
					break;
				default:
					throw std::invalid_argument(  "Do not recognize type of plaquette. Should be 1,2 or 3 (normal plaquette, temporal plaquette, spatial plaquette)" );
					break;
			}
		}
	private:
		int typeOfPlaquette;
	};

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

	class PolyakovloopTester : public GaugefieldTester {
	public:
		PolyakovloopTester(std::string inputfile):
			GaugefieldTester("polyakov", inputfile, 2) {
			const hardware::buffers::Plain<hmc_complex> pol(1, device);
			code->polyakov_device(getGaugefieldBuffer(), &pol);

			hmc_complex kernelResult_tmp;
			pol.dump(&kernelResult_tmp);
			kernelResult[0] = kernelResult_tmp.re;
			kernelResult[1] = kernelResult_tmp.im;
		}
	};

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

	class StoutSmearTester : public GaugefieldTester {
	public:
		StoutSmearTester(std::string inputfile):
			GaugefieldTester("stout_smear", inputfile) {
			auto gaugefieldBuffer = getGaugefieldBuffer();

			const hardware::buffers::Plain<hmc_float> plaq(1, device );
			const hardware::buffers::Plain<hmc_float> splaq(1, device);
			const hardware::buffers::Plain<hmc_float> tplaq(1, device);
			const hardware::buffers::SU3 out(gaugefieldBuffer->get_elements(), device);

			code->stout_smear_device( gaugefieldBuffer, &out);

			code->plaquette_device( &out, &plaq, &tplaq, &splaq);
			plaq.dump(&kernelResult[0]);
		}
	};

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

	class RectanglesTester : public GaugefieldTester {
	public:
		RectanglesTester(std::string inputfile):
			GaugefieldTester("rectangles", inputfile) {
			const hardware::buffers::Plain<hmc_float> rect(1, device );
			code->rectangles_device(getGaugefieldBuffer(), &rect);
			rect.dump(&kernelResult[0]);
		}
	};

	BOOST_AUTO_TEST_CASE( RECTANGLES_1 )
	{
		RectanglesTester rectanglesTester("rectangles_input_1");
	}

	BOOST_AUTO_TEST_CASE( RECTANGLES_2 )
	{
		RectanglesTester rectanglesTester("rectangles_input_2");
	}

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

