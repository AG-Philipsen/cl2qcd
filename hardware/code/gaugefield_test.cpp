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

#include "testUtilities.hpp"
#include "kernelTester.hpp"
#include "gaugefield.hpp"
#include "../../host_functionality/host_operations_gaugefield.h"
#include "mockups.hpp"

const int ns4 = 4;
const int nt4 = 4;
const int ns8 = 8;
const int nt8 = 8;
const int ns12 = 12;
const int nt12 = 12;
const int ns16 = 16;
const int nt16 = 16;

enum FillType {cold = 1, nonTrivial};

struct TestParameters {
	std::vector<double> referenceValue;
	size_t numberOfValues;
	int typeOfComparison;

	TestParameters(std::vector<double> referenceValueIn){
		referenceValue = referenceValueIn;
		numberOfValues = referenceValueIn.size();
		typeOfComparison = 1;
	}

	TestParameters(std::vector<double> referenceValueIn, int typeOfComparisonIn){
		referenceValue = referenceValueIn;
		numberOfValues = referenceValueIn.size();
		typeOfComparison = typeOfComparisonIn;
	}
};

//todo: check if this can be moved elsewhere...
void set_cold(Matrixsu3 * field, size_t elems)
{
	for(size_t i = 0; i < elems; ++i) {
		field[i] = unit_matrixsu3();
	}
}

void setGaugefield(Matrixsu3 * field, size_t elems, const FillType fillTypeIn)
{
	for(size_t i = 0; i < elems; ++i)
	{
		switch (fillTypeIn)
		{
		case cold:
			field[i] = unit_matrixsu3();
			break;
		case nonTrivial:
			field[i] = nonTrivialSu3Matrix();
			break;
		default:
			BOOST_ERROR("No valid FillType specified");
		}
	}
}

const Matrixsu3* createGaugefield(const int numberOfElements, const FillType fillTypeIn)
{
	Matrixsu3 * tmp = new Matrixsu3[numberOfElements];
	setGaugefield( tmp, numberOfElements, fillTypeIn);
	return tmp;
}

class GaugefieldTester : public KernelTester {
public:
		GaugefieldTester(std::string kernelName, const hardware::HardwareParametersMockup & hardwareParameters, const hardware::code::OpenClKernelParametersMockup & kernelParameters,
				const hardware::OpenClCodeMockup & kernelBuilder, struct TestParameters testParams, const FillType fillTypeIn):
			KernelTester(kernelName, hardwareParameters, kernelParameters, kernelBuilder, testParams.numberOfValues, testParams.typeOfComparison, testParams.referenceValue ),
			numberOfElements(hardwareParameters.getLatticeVolume() * NDIM)
	{
		gaugefieldBuffer = new hardware::buffers::SU3( numberOfElements, device);
		const Matrixsu3 * gf_host = createGaugefield(numberOfElements, fillTypeIn);
		device->getGaugefieldCode()->importGaugefield(gaugefieldBuffer, gf_host);
		delete[] gf_host;


		code = device->getGaugefieldCode();
	}

	~GaugefieldTester()
	{
		delete gaugefieldBuffer;
	}

protected:
	const int numberOfElements;
	const hardware::code::Gaugefield * code;
	const hardware::buffers::SU3 * gaugefieldBuffer;

};

BOOST_AUTO_TEST_SUITE ( PLAQUETTE )

	class PlaquetteTester : public GaugefieldTester {
	public:
			PlaquetteTester(const hardware::HardwareParametersMockup & hardwareParameters, const hardware::code::OpenClKernelParametersMockup & kernelParameters,
							const hardware::OpenClCodeMockup & kernelBuilder, struct TestParameters testParams, const FillType fillTypeIn, int typeOfPlaquette):
				GaugefieldTester("plaquette", hardwareParameters, kernelParameters, kernelBuilder, testParams, fillTypeIn), typeOfPlaquette(typeOfPlaquette)  {
			const hardware::buffers::Plain<hmc_float> plaq(1, device);
			const hardware::buffers::Plain<hmc_float> splaq(1, device);
			const hardware::buffers::Plain<hmc_float> tplaq(1, device);
			code->plaquette_device(gaugefieldBuffer, &plaq, &tplaq, &splaq);

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
		std::vector<double> referenceValue = {1536.};
		hardware::HardwareParametersMockup hardwareParameters(ns4,nt4);
		hardware::code::OpenClKernelParametersMockup kernelParameters(ns4,nt4);
		hardware::OpenClCodeMockup kernelBuilder(kernelParameters);
		TestParameters testParameters {referenceValue};
		PlaquetteTester plaquetteTester(hardwareParameters, kernelParameters, kernelBuilder, testParameters, FillType::cold, 1);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_2 )
	{
		std::vector<double> referenceValue = {1536.002605};
			hardware::HardwareParametersMockup hardwareParameters(ns4,nt4);
			hardware::code::OpenClKernelParametersMockup kernelParameters(ns4,nt4);
			hardware::OpenClCodeMockup kernelBuilder(kernelParameters);
		TestParameters testParameters {referenceValue};
		PlaquetteTester plaquetteTester(hardwareParameters, kernelParameters, kernelBuilder, testParameters, FillType::nonTrivial, 1);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_TEMPORAL_1 )
	{
		std::vector<double> referenceValue = {768.};
		hardware::HardwareParametersMockup hardwareParameters(ns4,nt4);
		hardware::code::OpenClKernelParametersMockup kernelParameters(ns4,nt4);
		hardware::OpenClCodeMockup kernelBuilder(kernelParameters);
		TestParameters testParameters {referenceValue};
		PlaquetteTester plaquetteTester(hardwareParameters, kernelParameters, kernelBuilder, testParameters, FillType::cold,  2);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_TEMPORAL_2 )
	{
		std::vector<double> referenceValue = {768.00130250240136};
		hardware::HardwareParametersMockup hardwareParameters(ns4,nt4);
		hardware::code::OpenClKernelParametersMockup kernelParameters(ns4,nt4);
		hardware::OpenClCodeMockup kernelBuilder(kernelParameters);
		TestParameters testParameters {referenceValue};
		PlaquetteTester plaquetteTester(hardwareParameters, kernelParameters, kernelBuilder, testParameters, FillType::nonTrivial, 2);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_SPATIAL_1 )
	{
		std::vector<double> referenceValue = {768.};
		hardware::HardwareParametersMockup hardwareParameters(ns4,nt4);
		hardware::code::OpenClKernelParametersMockup kernelParameters(ns4,nt4);
		hardware::OpenClCodeMockup kernelBuilder(kernelParameters);
		TestParameters testParameters {referenceValue};
		PlaquetteTester plaquetteTester(hardwareParameters, kernelParameters, kernelBuilder, testParameters, FillType::cold, 3);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_SPATIAL_2 )
	{
		std::vector<double> referenceValue = {768.00130250240136};
		hardware::HardwareParametersMockup hardwareParameters(ns4,nt4);
		hardware::code::OpenClKernelParametersMockup kernelParameters(ns4,nt4);
		hardware::OpenClCodeMockup kernelBuilder(kernelParameters);
		TestParameters testParameters {referenceValue};
		PlaquetteTester plaquetteTester(hardwareParameters, kernelParameters, kernelBuilder, testParameters, FillType::nonTrivial,  3);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_REDUCTION_1 )
	{
		std::vector<double> referenceValue = {24576};
		hardware::HardwareParametersMockup hardwareParameters(ns8,nt8);
		hardware::code::OpenClKernelParametersMockup kernelParameters(ns8,nt8);
		hardware::OpenClCodeMockup kernelBuilder(kernelParameters);
		TestParameters testParameters {referenceValue};
		PlaquetteTester plaquetteTester(hardwareParameters, kernelParameters, kernelBuilder, testParameters, FillType::cold,  1);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_REDUCTION_2 )
	{
		std::vector<double> referenceValue = {124416};
		hardware::HardwareParametersMockup hardwareParameters(ns12,nt12);
		hardware::code::OpenClKernelParametersMockup kernelParameters(ns12,nt12);
		hardware::OpenClCodeMockup kernelBuilder(kernelParameters);
		TestParameters testParameters {referenceValue};
		PlaquetteTester plaquetteTester(hardwareParameters, kernelParameters, kernelBuilder, testParameters, FillType::cold,  1);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_REDUCTION_3 )
	{
		std::vector<double> referenceValue = {393216};
		hardware::HardwareParametersMockup hardwareParameters(ns16,nt16);
		hardware::code::OpenClKernelParametersMockup kernelParameters(ns16,nt16);
		hardware::OpenClCodeMockup kernelBuilder(kernelParameters);
		TestParameters testParameters {referenceValue};
		PlaquetteTester plaquetteTester(hardwareParameters, kernelParameters, kernelBuilder, testParameters, FillType::cold,  1);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( POLYAKOV )

	class PolyakovloopTester : public GaugefieldTester {
	public:
			PolyakovloopTester(const hardware::HardwareParametersMockup & hardwareParameters, const hardware::code::OpenClKernelParametersMockup & kernelParameters,
					const hardware::OpenClCodeMockup & kernelBuilder, struct TestParameters testParams, const FillType fillTypeIn):
				GaugefieldTester("polyakov", hardwareParameters, kernelParameters, kernelBuilder, testParams, fillTypeIn) {
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
		std::vector<double> referenceValue = {64.,0.};
		hardware::HardwareParametersMockup hardwareParameters(ns4,nt4);
		hardware::code::OpenClKernelParametersMockup kernelParameters(ns4,nt4);
		hardware::OpenClCodeMockup kernelBuilder(kernelParameters);
		TestParameters testParameters {referenceValue};
		PolyakovloopTester polyakovloopTester(hardwareParameters, kernelParameters, kernelBuilder, testParameters, FillType::cold);
	}

	BOOST_AUTO_TEST_CASE( POLYAKOV_2)
	{
		std::vector<double> referenceValue = {-17.1117721375,-31.0747993518};
		hardware::HardwareParametersMockup hardwareParameters(ns4,nt4);
		hardware::code::OpenClKernelParametersMockup kernelParameters(ns4,nt4);
		hardware::OpenClCodeMockup kernelBuilder(kernelParameters);
		TestParameters testParameters {referenceValue};
		PolyakovloopTester polyakovloopTester(hardwareParameters, kernelParameters, kernelBuilder, testParameters, FillType::nonTrivial);
	}

	BOOST_AUTO_TEST_CASE( POLYAKOV_REDUCTION_1 )
	{
		std::vector<double> referenceValue = {512.,0.};
		hardware::HardwareParametersMockup hardwareParameters(ns8,nt8);
		hardware::code::OpenClKernelParametersMockup kernelParameters(ns8,nt8);
		hardware::OpenClCodeMockup kernelBuilder(kernelParameters);
		TestParameters testParameters {referenceValue};
		PolyakovloopTester polyakovloopTester(hardwareParameters, kernelParameters, kernelBuilder, testParameters, FillType::cold);
	}

	BOOST_AUTO_TEST_CASE( POLYAKOV_REDUCTION_2 )
	{
		std::vector<double> referenceValue = {1728.,0.};
		hardware::HardwareParametersMockup hardwareParameters(ns12,nt12);
		hardware::code::OpenClKernelParametersMockup kernelParameters(ns12,nt12);
		hardware::OpenClCodeMockup kernelBuilder(kernelParameters);
		TestParameters testParameters {referenceValue};
		PolyakovloopTester polyakovloopTester(hardwareParameters, kernelParameters, kernelBuilder, testParameters, FillType::cold);
	}

	BOOST_AUTO_TEST_CASE( POLYAKOV_REDUCTION_3 )
	{
		std::vector<double> referenceValue = {4096.,0.};
		hardware::HardwareParametersMockup hardwareParameters(ns16,nt16);
		hardware::code::OpenClKernelParametersMockup kernelParameters(ns16,nt16);
		hardware::OpenClCodeMockup kernelBuilder(kernelParameters);
		TestParameters testParameters {referenceValue};
		PolyakovloopTester polyakovloopTester(hardwareParameters, kernelParameters, kernelBuilder, testParameters, FillType::cold);
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
			StoutSmearTester(const hardware::HardwareParametersMockup & hardwareParameters, const hardware::code::OpenClKernelParametersMockup & kernelParameters,
					const hardware::OpenClCodeMockup & kernelBuilder, struct TestParameters testParams, const FillType fillTypeIn):
				GaugefieldTester("stout_smear", hardwareParameters, kernelParameters, kernelBuilder, testParams, fillTypeIn) {

			const hardware::buffers::Plain<hmc_float> plaq(1, device );
			const hardware::buffers::Plain<hmc_float> splaq(1, device);
			const hardware::buffers::Plain<hmc_float> tplaq(1, device);
			const hardware::buffers::SU3 out(gaugefieldBuffer->get_elements(), device);

			code->stout_smear_device( gaugefieldBuffer, &out);

			code->plaquette_device( &out, &plaq, &tplaq, &splaq);
			plaq.dump(&kernelResult[0]);
		}
	};

//	BOOST_AUTO_TEST_CASE( STOUT_SMEAR_1 )
//	{
//		hardware::HardwareParametersMockup hardwareParameters(4,4);
//		FillType fillType = FillType::cold;
//		std::vector<std::string> parameterStrings {"--startcondition=continue", "--sourcefile=conf.00200", "--use_smearing=true", "--rho_iter=1", "--rho=0.001", "--prec=64", "--nspace=4", "--ntime=4", "--fermact=TWISTEDMASS", "--gaugeact=wilson", "--beta=5.69", "--solver_prec=1e-8",  "--test_ref_val=882.11113688812929"};
//		StoutSmearTester StoutSmearTester(parameterStrings, hardwareParameters, fillType);
//	}
//
////	BOOST_AUTO_TEST_CASE( STOUT_SMEAR_2 )
////	{
////		hardware::HardwareParametersMockup hardwareParameters(4,4);
////		FillType fillType = FillType::cold;
////		std::vector<std::string> parameterStrings {"--startcondition=continue", "--sourcefile=conf.00200", "--use_smearing=true", "--rho_iter=1", "--rho=0.", "--prec=64", "--nspace=4", "--ntime=4", "--fermact=TWISTEDMASS", "--gaugeact=wilson", "--beta=5.69", "--solver_prec=1e-8",  "--test_ref_val=877.17444356279361"};
////		StoutSmearTester StoutSmearTester(parameterStrings, hardwareParameters, fillType);
////	}
////
////	BOOST_AUTO_TEST_CASE( STOUT_SMEAR_3 )
////	{
////		hardware::HardwareParametersMockup hardwareParameters(4,4);
////		FillType fillType = FillType::cold;
////		std::vector<std::string> parameterStrings {"--startcondition=continue", "--sourcefile=conf.00200", "--use_smearing=true", "--rho_iter=1", "--rho=0.001538", "--prec=64", "--nspace=4", "--ntime=4", "--fermact=TWISTEDMASS", "--gaugeact=wilson", "--beta=5.69", "--solver_prec=1e-8",  "--test_ref_val=884.76195718059716"};
////		StoutSmearTester StoutSmearTester(parameterStrings, hardwareParameters, fillType);
////	}
////
	BOOST_AUTO_TEST_CASE( STOUT_SMEAR_4 )
	{
		const int rhoIterIn = 1;
		const double rhoIn = 0.;
		const bool useSmearingIn = true;
		std::vector<double> referenceValue = {1536};
		hardware::HardwareParametersMockup hardwareParameters(ns4,nt4);
		hardware::code::OpenClKernelParametersMockup kernelParameters(ns4,nt4,rhoIterIn,rhoIn,useSmearingIn);
		hardware::OpenClCodeMockup kernelBuilder(kernelParameters);
		TestParameters testParameters {referenceValue};
		StoutSmearTester StoutSmearTester(hardwareParameters, kernelParameters, kernelBuilder, testParameters, FillType::cold);
	}

	BOOST_AUTO_TEST_CASE( STOUT_SMEAR_5 )
	{
		const int rhoIterIn = 1;
		const double rhoIn = 0.001;
		const bool useSmearingIn = true;
		std::vector<double> referenceValue = {1536};
		hardware::HardwareParametersMockup hardwareParameters(ns4,nt4);
		hardware::code::OpenClKernelParametersMockup kernelParameters(ns4,nt4,rhoIterIn,rhoIn,useSmearingIn);
		hardware::OpenClCodeMockup kernelBuilder(kernelParameters);
		TestParameters testParameters {referenceValue};
		StoutSmearTester StoutSmearTester(hardwareParameters, kernelParameters, kernelBuilder, testParameters, FillType::cold);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( RECTANGLES )

	class RectanglesTester : public GaugefieldTester {
	public:
			RectanglesTester(const hardware::HardwareParametersMockup & hardwareParameters, const hardware::code::OpenClKernelParametersMockup & kernelParameters,
					const hardware::OpenClCodeMockup & kernelBuilder, struct TestParameters testParams, const FillType fillTypeIn):
				GaugefieldTester("rectangles", hardwareParameters, kernelParameters, kernelBuilder, testParams, fillTypeIn) {
			const hardware::buffers::Plain<hmc_float> rect(1, device );
			code->rectangles_device(gaugefieldBuffer, &rect);
			rect.dump(&kernelResult[0]);
		}
	};

	BOOST_AUTO_TEST_CASE( RECTANGLES_1 )
	{
		const bool useRectanglesIn = true;
		std::vector<double> referenceValue = {3072.};
		hardware::HardwareParametersMockup hardwareParameters(ns4,nt4);
		hardware::code::OpenClKernelParametersMockup kernelParameters(ns4,nt4,useRectanglesIn);
		hardware::OpenClCodeMockup kernelBuilder(kernelParameters);
		TestParameters testParameters {referenceValue};
		RectanglesTester rectanglesTester(hardwareParameters, kernelParameters, kernelBuilder, testParameters,FillType::cold);
	}

	BOOST_AUTO_TEST_CASE( RECTANGLES_2 )
	{
		const bool useRectanglesIn = true;
		std::vector<double> referenceValue = {3072.00781501};
		hardware::HardwareParametersMockup hardwareParameters(ns4,nt4);
		hardware::code::OpenClKernelParametersMockup kernelParameters(ns4,nt4,useRectanglesIn);
		hardware::OpenClCodeMockup kernelBuilder(kernelParameters);
		TestParameters testParameters {referenceValue};
		RectanglesTester rectanglesTester(hardwareParameters, kernelParameters, kernelBuilder, testParameters,FillType::nonTrivial);
	}

	BOOST_AUTO_TEST_CASE( RECTANGLES_REDUCTION_1 )
	{
			const bool useRectanglesIn = true;
			std::vector<double> referenceValue = {49152.};
			hardware::HardwareParametersMockup hardwareParameters(ns8,nt8);
			hardware::code::OpenClKernelParametersMockup kernelParameters(ns8,nt8,useRectanglesIn);
			hardware::OpenClCodeMockup kernelBuilder(kernelParameters);
			TestParameters testParameters {referenceValue};
			RectanglesTester rectanglesTester(hardwareParameters, kernelParameters, kernelBuilder, testParameters,FillType::cold);
	}

	BOOST_AUTO_TEST_CASE( RECTANGLES_REDUCTION_2 )
	{
		const bool useRectanglesIn = true;
		std::vector<double> referenceValue = {248832.};
		hardware::HardwareParametersMockup hardwareParameters(ns12,nt12);
		hardware::code::OpenClKernelParametersMockup kernelParameters(ns12,nt12,useRectanglesIn);
		hardware::OpenClCodeMockup kernelBuilder(kernelParameters);
		TestParameters testParameters {referenceValue};
		RectanglesTester rectanglesTester(hardwareParameters, kernelParameters, kernelBuilder, testParameters,FillType::cold);
	}

	BOOST_AUTO_TEST_CASE( RECTANGLES_REDUCTION_3 )
	{
		const bool useRectanglesIn = true;
		std::vector<double> referenceValue = {786432.};
		hardware::HardwareParametersMockup hardwareParameters(ns16,nt16);
		hardware::code::OpenClKernelParametersMockup kernelParameters(ns16,nt16,useRectanglesIn);
		hardware::OpenClCodeMockup kernelBuilder(kernelParameters);
		TestParameters testParameters {referenceValue};
		RectanglesTester rectanglesTester(hardwareParameters, kernelParameters, kernelBuilder, testParameters,FillType::cold);
	}

BOOST_AUTO_TEST_SUITE_END()

