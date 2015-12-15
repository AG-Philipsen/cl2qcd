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

template<typename TesterClass, typename ParameterClass> void callTest(const ParameterClass parametersForThisTest)
{
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns,parametersForThisTest.nt);
	hardware::code::OpenClKernelParametersMockup kernelParameters(parametersForThisTest.ns,parametersForThisTest.nt);
	ParameterCollection parameterCollection(hardwareParameters, kernelParameters);
	TesterClass tester(parameterCollection, parametersForThisTest);
}

template<typename ParameterClass, typename TesterClass> void performTest(const ReferenceValues refValuesIn, const LatticeExtents latticeExtentsIn, const GaugefieldFillType fillType)
{
	ParameterClass testParams {refValuesIn, latticeExtentsIn, fillType};
	callTest<TesterClass, ParameterClass>(testParams);
}

template<typename ParameterClass, typename TesterClass, typename additionalArgument> void performTest(const ReferenceValues refValuesIn, const LatticeExtents latticeExtentsIn, const GaugefieldFillType fillType, const additionalArgument addArg)
{
	ParameterClass testParams {refValuesIn, latticeExtentsIn, fillType, addArg};
	callTest<TesterClass, ParameterClass>(testParams);
}

template<typename ParameterClass, typename TesterClass> void performTest(const ReferenceValues refValuesIn, const LatticeExtents latticeExtentsIn, const GaugefieldFillType fillType, const int rhoIter, const double rho)
{
	const bool useSmearing = true;
	ParameterClass testParams {refValuesIn, latticeExtentsIn, fillType, rhoIter, rho};
	hardware::HardwareParametersMockup hardwareParameters(testParams.ns,testParams.nt);
	hardware::code::OpenClKernelParametersMockup kernelParameters(testParams.ns,testParams.nt, testParams.rhoIter, testParams.rho, useSmearing);
	ParameterCollection parameterCollection(hardwareParameters, kernelParameters);
	TesterClass tester(parameterCollection, testParams);
}

template<typename TesterClass> void performTest(const ReferenceValues refValuesIn, const LatticeExtents latticeExtentsIn, const GaugefieldFillType fillType)
{
	const bool useRectangles = true;
	GaugefieldTestParameters testParams {refValuesIn, latticeExtentsIn, fillType};
	hardware::HardwareParametersMockup hardwareParameters(testParams.ns,testParams.nt);
	hardware::code::OpenClKernelParametersMockup kernelParameters(testParams.ns,testParams.nt, useRectangles);
	ParameterCollection parameterCollection(hardwareParameters, kernelParameters);
	TesterClass tester(parameterCollection, testParams);
}

BOOST_AUTO_TEST_SUITE ( PLAQUETTE )

	enum TypeOfPlaquette {plaquette = 1, temporalPlaquette, spatialPlaquette };

	struct PlaquetteTestParameters : public GaugefieldTestParameters
	{
		PlaquetteTestParameters(const ReferenceValues referenceValueIn, const LatticeExtents latticeExtentsIn, const GaugefieldFillType fillTypeIn, const TypeOfPlaquette typeOfPlaquetteIn) :
			GaugefieldTestParameters(referenceValueIn,latticeExtentsIn, fillTypeIn), typeOfPlaquette(typeOfPlaquetteIn), TestParameters(latticeExtentsIn){}
		TypeOfPlaquette typeOfPlaquette;
	};

	class PlaquetteTester : public GaugefieldTester {
	public:
		PlaquetteTester(const ParameterCollection & parameterCollection, const PlaquetteTestParameters testParams):
			GaugefieldTester("plaquette", parameterCollection,  testParams), typeOfPlaquette(testParams.typeOfPlaquette)
		{
			const hardware::buffers::Plain<hmc_float> plaq(1, device);
			const hardware::buffers::Plain<hmc_float> splaq(1, device);
			const hardware::buffers::Plain<hmc_float> tplaq(1, device);
			code->plaquette_device(gaugefieldBuffer, &plaq, &tplaq, &splaq);

			switch( typeOfPlaquette ) {
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
	private:
		TypeOfPlaquette typeOfPlaquette;
	};

	BOOST_AUTO_TEST_CASE( PLAQUETTE_1 )
	{
		performTest<PlaquetteTestParameters, PlaquetteTester, TypeOfPlaquette>(ReferenceValues{1536.}, LatticeExtents{ns4, nt4}, GaugefieldFillType::cold, TypeOfPlaquette::plaquette );
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_2 )
	{
		performTest<PlaquetteTestParameters, PlaquetteTester, TypeOfPlaquette>(ReferenceValues{1536.002605}, LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial, TypeOfPlaquette::plaquette );
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_TEMPORAL_1 )
	{
		performTest<PlaquetteTestParameters, PlaquetteTester, TypeOfPlaquette>(ReferenceValues{768.}, LatticeExtents{ns4, nt4}, GaugefieldFillType::cold, TypeOfPlaquette::temporalPlaquette);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_TEMPORAL_2 )
	{
		performTest<PlaquetteTestParameters, PlaquetteTester, TypeOfPlaquette>(ReferenceValues{768.00130250240136}, LatticeExtents{LatticeExtents{ns4, nt4}}, GaugefieldFillType::nonTrivial, TypeOfPlaquette::temporalPlaquette);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_SPATIAL_1 )
	{
		performTest<PlaquetteTestParameters, PlaquetteTester, TypeOfPlaquette>(ReferenceValues{768.}, LatticeExtents{LatticeExtents{ns4, nt4}}, GaugefieldFillType::cold, TypeOfPlaquette::temporalPlaquette);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_SPATIAL_2 )
	{
		performTest<PlaquetteTestParameters, PlaquetteTester, TypeOfPlaquette>(ReferenceValues{768.00130250240136}, LatticeExtents{LatticeExtents{ns4, nt4}}, GaugefieldFillType::nonTrivial, TypeOfPlaquette::spatialPlaquette);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_REDUCTION_1 )
	{
		performTest<PlaquetteTestParameters, PlaquetteTester, TypeOfPlaquette>(ReferenceValues{24576}, LatticeExtents{ns8, nt8}, GaugefieldFillType::cold, TypeOfPlaquette::plaquette);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_REDUCTION_2 )
	{
		performTest<PlaquetteTestParameters, PlaquetteTester, TypeOfPlaquette>(ReferenceValues{124416}, LatticeExtents{ns12, nt12}, GaugefieldFillType::cold, TypeOfPlaquette::plaquette);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_REDUCTION_3 )
	{
		performTest<PlaquetteTestParameters, PlaquetteTester, TypeOfPlaquette>(ReferenceValues{393216}, LatticeExtents{ns16, nt16}, GaugefieldFillType::cold, TypeOfPlaquette::plaquette);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( POLYAKOV )

	class PolyakovloopTester : public GaugefieldTester {
	public:
			PolyakovloopTester(const ParameterCollection & parameterCollection,
					GaugefieldTestParameters testParams):
				GaugefieldTester("polyakov", parameterCollection,  testParams) {
			const hardware::buffers::Plain<hmc_complex> pol(1, device);
			code->polyakov_device(gaugefieldBuffer, &pol);

			hmc_complex kernelResult_tmp;
			pol.dump(&kernelResult_tmp);
			kernelResult[0] = kernelResult_tmp.re;
			kernelResult[1] = kernelResult_tmp.im;
		}
	};

	BOOST_AUTO_TEST_CASE( POLYAKOV_1 )
	{
		performTest<GaugefieldTestParameters, PolyakovloopTester>(ReferenceValues{64.,0.}, LatticeExtents{LatticeExtents{ns4, nt4}}, GaugefieldFillType::cold);
	}

	BOOST_AUTO_TEST_CASE( POLYAKOV_2)
	{
		performTest<GaugefieldTestParameters, PolyakovloopTester>(ReferenceValues{-17.1117721375,-31.0747993518}, LatticeExtents{LatticeExtents{ns4, nt4}}, GaugefieldFillType::nonTrivial);
	}

	BOOST_AUTO_TEST_CASE( POLYAKOV_REDUCTION_1 )
	{
		performTest<GaugefieldTestParameters, PolyakovloopTester>(ReferenceValues{512.,0.}, LatticeExtents{ns8, nt8}, GaugefieldFillType::cold);
	}

	BOOST_AUTO_TEST_CASE( POLYAKOV_REDUCTION_2 )
	{
		performTest<GaugefieldTestParameters, PolyakovloopTester>(ReferenceValues{1728.,0.}, LatticeExtents{ns12, nt12}, GaugefieldFillType::cold);
	}

	BOOST_AUTO_TEST_CASE( POLYAKOV_REDUCTION_3 )
	{
		performTest<GaugefieldTestParameters, PolyakovloopTester>(ReferenceValues{4096.,0.}, LatticeExtents{ns16, nt16}, GaugefieldFillType::cold);
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

	struct SmearingTestParameters : public GaugefieldTestParameters
	{
		int rhoIter;
		double rho;

		SmearingTestParameters(std::vector<double> referenceValueIn, const LatticeExtents latticeExtentsIn, GaugefieldFillType fillTypeIn, int rhoIterIn, double rhoIn) :
			GaugefieldTestParameters(referenceValueIn,latticeExtentsIn, fillTypeIn), rhoIter(rhoIterIn), rho(rhoIn), TestParameters(latticeExtentsIn) {}
	};

	class StoutSmearTester : public GaugefieldTester {
	public:
			StoutSmearTester(const ParameterCollection & parameterCollection,
					struct SmearingTestParameters testParams):
				GaugefieldTester("stout_smear", parameterCollection,  testParams) {

			const hardware::buffers::Plain<hmc_float> plaq(1, device );
			const hardware::buffers::Plain<hmc_float> splaq(1, device);
			const hardware::buffers::Plain<hmc_float> tplaq(1, device);
			const hardware::buffers::SU3 out(gaugefieldBuffer->get_elements(), device);

			code->stout_smear_device( gaugefieldBuffer, &out);

			code->plaquette_device( &out, &plaq, &tplaq, &splaq);
			plaq.dump(&kernelResult[0]);
		}
	};

	BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES (STOUT_SMEAR_1, 1)
	BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES (STOUT_SMEAR_2, 1)
	BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES (STOUT_SMEAR_3, 1)

	BOOST_AUTO_TEST_CASE( STOUT_SMEAR_1 )
	{
		const int rhoIterIn = 1;
		const double rhoIn = 0.001;
		performTest<SmearingTestParameters, StoutSmearTester>( ReferenceValues{-1234}, LatticeExtents{LatticeExtents{ns4, nt4}}, GaugefieldFillType::nonTrivial, rhoIterIn, rhoIn);
	}

	BOOST_AUTO_TEST_CASE( STOUT_SMEAR_2 )
	{
		const int rhoIterIn = 1;
		const double rhoIn = 0.;
		performTest<SmearingTestParameters, StoutSmearTester>( ReferenceValues{-1234}, LatticeExtents{LatticeExtents{ns4, nt4}}, GaugefieldFillType::nonTrivial, rhoIterIn, rhoIn);
	}

	BOOST_AUTO_TEST_CASE( STOUT_SMEAR_3 )
	{
		const int rhoIterIn = 1;
		const double rhoIn = 0.001538;
		performTest<SmearingTestParameters, StoutSmearTester>( ReferenceValues{-1234}, LatticeExtents{LatticeExtents{ns4, nt4}}, GaugefieldFillType::nonTrivial, rhoIterIn, rhoIn);
	}

	BOOST_AUTO_TEST_CASE( STOUT_SMEAR_4 )
	{
		const int rhoIterIn = 1;
		const double rhoIn = 0.;
		performTest<SmearingTestParameters, StoutSmearTester>( ReferenceValues{1536}, LatticeExtents{LatticeExtents{ns4, nt4}}, GaugefieldFillType::cold, rhoIterIn, rhoIn);
	}

	BOOST_AUTO_TEST_CASE( STOUT_SMEAR_5 )
	{
		const int rhoIterIn = 1;
		const double rhoIn = 0.001;
		performTest<SmearingTestParameters, StoutSmearTester>( ReferenceValues{1536}, LatticeExtents{LatticeExtents{ns4, nt4}}, GaugefieldFillType::cold, rhoIterIn, rhoIn);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( RECTANGLES )

	class RectanglesTester : public GaugefieldTester {
	public:
		RectanglesTester(const ParameterCollection & parameterCollection, const GaugefieldTestParameters testParams):
			GaugefieldTester("rectangles", parameterCollection,  testParams)
		{
			const hardware::buffers::Plain<hmc_float> rect(1, device );
			code->rectangles_device(gaugefieldBuffer, &rect);
			rect.dump(&kernelResult[0]);
		}
	};

	BOOST_AUTO_TEST_CASE( RECTANGLES_1 )
	{
		performTest<RectanglesTester>(ReferenceValues{3072.}, LatticeExtents{ns4, nt4}, GaugefieldFillType::cold);
	}

	BOOST_AUTO_TEST_CASE( RECTANGLES_2 )
	{
		performTest<RectanglesTester>(ReferenceValues{3072.00781501}, LatticeExtents{ns4, nt4}, GaugefieldFillType::nonTrivial);
	}

	BOOST_AUTO_TEST_CASE( RECTANGLES_REDUCTION_1 )
	{
		performTest<RectanglesTester>(ReferenceValues{49152.}, LatticeExtents{ns8, nt8}, GaugefieldFillType::cold);
	}

	BOOST_AUTO_TEST_CASE( RECTANGLES_REDUCTION_2 )
	{
		performTest<RectanglesTester>(ReferenceValues{248832.}, LatticeExtents{ns12, nt12}, GaugefieldFillType::cold);
	}

	BOOST_AUTO_TEST_CASE( RECTANGLES_REDUCTION_3 )
	{
		performTest<RectanglesTester>(ReferenceValues{786432.}, LatticeExtents{ns16, nt16}, GaugefieldFillType::cold);
	}

BOOST_AUTO_TEST_SUITE_END()

