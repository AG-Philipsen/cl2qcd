/*
 * Copyright 2012, 2013, 2014, 2015
 * 	Christopher Pinke, Matthias Bach, Francesca Cuteri
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

#include "GaugefieldTester.hpp"

enum TypeOfPlaquette {plaquette = 1, temporalPlaquette, spatialPlaquette };

struct SmearingTestParameters : public GaugefieldTestParameters
{
	int rhoIter;
	double rho;

	SmearingTestParameters(const LatticeExtents latticeExtentsIn, GaugefieldFillType fillTypeIn, int rhoIterIn, double rhoIn) :
		TestParameters(latticeExtentsIn), GaugefieldTestParameters(latticeExtentsIn, fillTypeIn), rhoIter(rhoIterIn), rho(rhoIn) {}
};

struct PlaquetteTester : public GaugefieldTester
{
	PlaquetteTester(const ParameterCollection & parameterCollection, const GaugefieldTestParameters testParams, const ReferenceValues rV, const TypeOfPlaquette typeOfPlaquette):
		GaugefieldTester("plaquette", parameterCollection,  testParams, rV)
	{
		const hardware::buffers::Plain<hmc_float> plaq(1, device);
		const hardware::buffers::Plain<hmc_float> splaq(1, device);
		const hardware::buffers::Plain<hmc_float> tplaq(1, device);
		GaugefieldCreator gf(testParams.latticeExtents);
		gaugefieldBuffer = new hardware::buffers::SU3(testParams.latticeExtents, this->device);
		const Matrixsu3 * gf_host = gf.createGaugefield(testParams.fillType);
        device->getGaugefieldCode()->importGaugefield(gaugefieldBuffer, gf_host);
        delete[] gf_host;

		code->plaquette_device(gaugefieldBuffer, &plaq, &tplaq, &splaq);

		switch( typeOfPlaquette )
		{
			case TypeOfPlaquette::plaquette:
				plaq.dump(&kernelResult[0]);
				break;
			case TypeOfPlaquette::temporalPlaquette:
				tplaq.dump(&kernelResult[0]);
				break;
			case TypeOfPlaquette::spatialPlaquette:
				splaq.dump(&kernelResult[0]);
				break;
			default:
				throw std::invalid_argument(  "Do not recognize type of plaquette" );
				break;
		}
	}
protected:
	const hardware::buffers::SU3 * gaugefieldBuffer;
};

struct PolyakovloopTester : public GaugefieldTester
{
	PolyakovloopTester(const ParameterCollection & parameterCollection, const GaugefieldTestParameters testParams, const ReferenceValues rV):
	GaugefieldTester("polyakov", parameterCollection,  testParams, rV)
	{
		const hardware::buffers::Plain<hmc_complex> pol(1, device);
		GaugefieldCreator gf(testParams.latticeExtents);
		gaugefieldBuffer = new hardware::buffers::SU3(testParams.latticeExtents, this->device);
		const Matrixsu3 * gf_host = gf.createGaugefield(testParams.fillType);
        device->getGaugefieldCode()->importGaugefield(gaugefieldBuffer, gf_host);
        delete[] gf_host;

		code->polyakov_device(gaugefieldBuffer, &pol);

		hmc_complex kernelResult_tmp;
		pol.dump(&kernelResult_tmp);
		kernelResult[0] = kernelResult_tmp.re;
		kernelResult[1] = kernelResult_tmp.im;
	}
protected:
	const hardware::buffers::SU3 * gaugefieldBuffer;
};

struct RectanglesTester : public GaugefieldTester
{
	RectanglesTester(const ParameterCollection & parameterCollection, const GaugefieldTestParameters testParams, const ReferenceValues rV):
		GaugefieldTester("rectangles", parameterCollection,  testParams, rV)
	{
		const hardware::buffers::Plain<hmc_float> rect(1, device );
		GaugefieldCreator gf(testParams.latticeExtents);
		gaugefieldBuffer = new hardware::buffers::SU3(testParams.latticeExtents, this->device);
		const Matrixsu3 * gf_host = gf.createGaugefield(testParams.fillType);
        device->getGaugefieldCode()->importGaugefield(gaugefieldBuffer, gf_host);
        delete[] gf_host;

		code->rectangles_device(gaugefieldBuffer, &rect);

		rect.dump(&kernelResult[0]);
	}
protected:
	const hardware::buffers::SU3 * gaugefieldBuffer;
};

struct StoutSmearTester : public GaugefieldTester
{
	StoutSmearTester(const ParameterCollection & parameterCollection, struct SmearingTestParameters testParams, const ReferenceValues rV):
		GaugefieldTester("stout_smear", parameterCollection,  testParams, rV)
	{
		const hardware::buffers::Plain<hmc_float> plaq(1, device );
		const hardware::buffers::Plain<hmc_float> splaq(1, device);
		const hardware::buffers::Plain<hmc_float> tplaq(1, device);
		GaugefieldCreator gf(testParams.latticeExtents);
		gaugefieldBuffer = new hardware::buffers::SU3(testParams.latticeExtents, this->device);
		const Matrixsu3 * gf_host = gf.createGaugefield(testParams.fillType);
        device->getGaugefieldCode()->importGaugefield(gaugefieldBuffer, gf_host);
        delete[] gf_host;
		const hardware::buffers::SU3 out(gaugefieldBuffer->get_elements(), device);

		code->stout_smear_device( gaugefieldBuffer, &out);
		code->plaquette_device( &out, &plaq, &tplaq, &splaq);

		plaq.dump(&kernelResult[0]);
	}
protected:
	const hardware::buffers::SU3 * gaugefieldBuffer;
};

template<typename TesterClass>
void performTest2(const ReferenceValues refValuesIn, const LatticeExtents latticeExtentsIn, const GaugefieldFillType fillType, const TypeOfPlaquette plaquetteType)
{
	GaugefieldTestParameters parametersForThisTest {latticeExtentsIn, fillType};
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns,parametersForThisTest.nt);
	hardware::code::OpenClKernelParametersMockup kernelParameters(parametersForThisTest.ns,parametersForThisTest.nt);
	ParameterCollection parameterCollection(hardwareParameters, kernelParameters);
	TesterClass tester(parameterCollection, parametersForThisTest, refValuesIn, plaquetteType);
}
template<typename TesterClass>
void performTest2(const ReferenceValues refValuesIn, const LatticeExtents latticeExtentsIn, const GaugefieldFillType fillType)
{
	GaugefieldTestParameters parametersForThisTest {latticeExtentsIn, fillType};
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns,parametersForThisTest.nt);
	hardware::code::OpenClKernelParametersMockup kernelParameters(parametersForThisTest.ns,parametersForThisTest.nt);
	ParameterCollection parameterCollection(hardwareParameters, kernelParameters);
	TesterClass tester(parameterCollection, parametersForThisTest, refValuesIn);
}

void performTest(const ReferenceValues refValuesIn, const LatticeExtents latticeExtentsIn, const GaugefieldFillType fillType, const int rhoIter, const double rho)
{
	const bool useSmearing = true;
	SmearingTestParameters testParams {latticeExtentsIn, fillType, rhoIter, rho};
	hardware::HardwareParametersMockup hardwareParameters(testParams.ns,testParams.nt);
	hardware::code::OpenClKernelParametersMockup kernelParameters(testParams.ns,testParams.nt, testParams.rhoIter, testParams.rho, useSmearing);
	ParameterCollection parameterCollection(hardwareParameters, kernelParameters);
	StoutSmearTester tester(parameterCollection, testParams, refValuesIn);
}

void testPlaquette(const ReferenceValues rV, const LatticeExtents lE, const GaugefieldFillType gF, const TypeOfPlaquette pT)
{
	performTest2<PlaquetteTester>(rV, lE, gF, pT );
}

void testPolyakov(const ReferenceValues rV, const LatticeExtents lE, const GaugefieldFillType gF)
{
	performTest2<PolyakovloopTester>(rV, lE, gF );
}

void testRectangles(const ReferenceValues rV, const LatticeExtents lE, const GaugefieldFillType gF)
{
	performTest2<RectanglesTester>(rV, lE, gF );
}

void testStoutSmear(const ReferenceValues rV, const LatticeExtents lE, const GaugefieldFillType gF, const double rhoIn)
{
	const int rhoIter = 1;
	performTest(rV, lE, gF, rhoIter, rhoIn );
}

BOOST_AUTO_TEST_SUITE ( PLAQUETTE )

	BOOST_AUTO_TEST_CASE( PLAQUETTE_1 )
	{
		testPlaquette(ReferenceValues{1536.}, LatticeExtents{ns4, nt4}, GaugefieldFillType::cold, TypeOfPlaquette::plaquette );
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_2 )
	{
		testPlaquette(ReferenceValues{1536.002605}, LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial, TypeOfPlaquette::plaquette );
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_TEMPORAL_1 )
	{
		testPlaquette(ReferenceValues{768.}, LatticeExtents{ns4, nt4}, GaugefieldFillType::cold, TypeOfPlaquette::temporalPlaquette);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_TEMPORAL_2 )
	{
		testPlaquette(ReferenceValues{768.00130250240136}, LatticeExtents{LatticeExtents{ns4, nt4}}, GaugefieldFillType::nonTrivial, TypeOfPlaquette::temporalPlaquette);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_SPATIAL_1 )
	{
		testPlaquette(ReferenceValues{768.}, LatticeExtents{LatticeExtents{ns4, nt4}}, GaugefieldFillType::cold, TypeOfPlaquette::temporalPlaquette);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_SPATIAL_2 )
	{
		testPlaquette(ReferenceValues{768.00130250240136}, LatticeExtents{LatticeExtents{ns4, nt4}}, GaugefieldFillType::nonTrivial, TypeOfPlaquette::spatialPlaquette);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_REDUCTION_1 )
	{
		testPlaquette(ReferenceValues{24576}, LatticeExtents{ns8, nt8}, GaugefieldFillType::cold, TypeOfPlaquette::plaquette);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_REDUCTION_2 )
	{
		testPlaquette(ReferenceValues{124416}, LatticeExtents{ns12, nt12}, GaugefieldFillType::cold, TypeOfPlaquette::plaquette);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_REDUCTION_3 )
	{
		testPlaquette(ReferenceValues{393216}, LatticeExtents{ns16, nt16}, GaugefieldFillType::cold, TypeOfPlaquette::plaquette);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( POLYAKOV )

	BOOST_AUTO_TEST_CASE( POLYAKOV_1 )
	{
		testPolyakov(ReferenceValues{64.,0.}, LatticeExtents{LatticeExtents{ns4, nt4}}, GaugefieldFillType::cold);
	}

	BOOST_AUTO_TEST_CASE( POLYAKOV_2)
	{
		testPolyakov(ReferenceValues{-17.1117721375,-31.0747993518}, LatticeExtents{LatticeExtents{ns4, nt4}}, GaugefieldFillType::nonTrivial);
	}

	BOOST_AUTO_TEST_CASE( POLYAKOV_REDUCTION_1 )
	{
		testPolyakov(ReferenceValues{512.,0.}, LatticeExtents{ns8, nt8}, GaugefieldFillType::cold);
	}

	BOOST_AUTO_TEST_CASE( POLYAKOV_REDUCTION_2 )
	{
		testPolyakov(ReferenceValues{1728.,0.}, LatticeExtents{ns12, nt12}, GaugefieldFillType::cold);
	}

	BOOST_AUTO_TEST_CASE( POLYAKOV_REDUCTION_3 )
	{
		testPolyakov(ReferenceValues{4096.,0.}, LatticeExtents{ns16, nt16}, GaugefieldFillType::cold);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( CONVERT_GAUGEFIELD_TO_SOA )

BOOST_AUTO_TEST_CASE( CONVERT_GAUGEFIELD_TO_SOA_1 )
{
	BOOST_TEST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( CONVERT_GAUGEFIELD_FROM_SOA )

BOOST_AUTO_TEST_CASE( CONVERT_GAUGEFIELD_FROM_SOA_1 )
{
	BOOST_TEST_MESSAGE("NOT YET IMPLEMENTED!");
}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( STOUT_SMEAR )

	BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES (STOUT_SMEAR_1, 1)
	BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES (STOUT_SMEAR_2, 1)
	BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES (STOUT_SMEAR_3, 1)

	BOOST_AUTO_TEST_CASE( STOUT_SMEAR_1 )
	{
		testStoutSmear( ReferenceValues{-1234}, LatticeExtents{LatticeExtents{ns4, nt4}}, GaugefieldFillType::nonTrivial, 0.001);
	}

	BOOST_AUTO_TEST_CASE( STOUT_SMEAR_2 )
	{
		testStoutSmear( ReferenceValues{-1234}, LatticeExtents{LatticeExtents{ns4, nt4}}, GaugefieldFillType::nonTrivial, 0.);
	}

	BOOST_AUTO_TEST_CASE( STOUT_SMEAR_3 )
	{
		testStoutSmear( ReferenceValues{-1234}, LatticeExtents{LatticeExtents{ns4, nt4}}, GaugefieldFillType::nonTrivial, 0.001538);
	}

	BOOST_AUTO_TEST_CASE( STOUT_SMEAR_4 )
	{
		testStoutSmear( ReferenceValues{1536}, LatticeExtents{LatticeExtents{ns4, nt4}}, GaugefieldFillType::cold, 0.);
	}

	BOOST_AUTO_TEST_CASE( STOUT_SMEAR_5 )
	{
		testStoutSmear( ReferenceValues{1536}, LatticeExtents{LatticeExtents{ns4, nt4}}, GaugefieldFillType::cold, 0.001);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( RECTANGLES )

	BOOST_AUTO_TEST_CASE( RECTANGLES_1 )
	{
		testRectangles(ReferenceValues{3072.}, LatticeExtents{ns4, nt4}, GaugefieldFillType::cold);
	}

	BOOST_AUTO_TEST_CASE( RECTANGLES_2 )
	{
		testRectangles(ReferenceValues{3072.00781501}, LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial);
	}

	BOOST_AUTO_TEST_CASE( RECTANGLES_REDUCTION_1 )
	{
		testRectangles(ReferenceValues{49152.}, LatticeExtents{ns8, nt8}, GaugefieldFillType::cold);
	}

	BOOST_AUTO_TEST_CASE( RECTANGLES_REDUCTION_2 )
	{
		testRectangles(ReferenceValues{248832.}, LatticeExtents{ns12, nt12}, GaugefieldFillType::cold);
	}

	BOOST_AUTO_TEST_CASE( RECTANGLES_REDUCTION_3 )
	{
		testRectangles(ReferenceValues{786432.}, LatticeExtents{ns16, nt16}, GaugefieldFillType::cold);
	}

BOOST_AUTO_TEST_SUITE_END()

