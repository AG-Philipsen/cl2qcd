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

struct TestParametersGaugeField : public TestParameters
{
	FillType fillType;

	TestParametersGaugeField(std::vector<double> referenceValueIn, int nsIn, int ntIn, FillType fillTypeIn):
		TestParameters(referenceValueIn, nsIn, ntIn) {
		referenceValue = referenceValueIn;
		ns = nsIn;
		nt = ntIn;
		fillType = fillTypeIn;
	}
};

struct TestParametersSpecToPlaq : public TestParametersGaugeField
{
	int typeOfPlaquette;

	TestParametersSpecToPlaq(std::vector<double> referenceValueIn, int nsIn, int ntIn, FillType fillTypeIn, int typeOfPlaquetteIn) :
		TestParametersGaugeField(referenceValueIn, nsIn, ntIn, fillTypeIn)
	{
		referenceValue = referenceValueIn;
		ns = nsIn;
		nt = ntIn;
		fillType = fillTypeIn;
		typeOfPlaquette = typeOfPlaquetteIn;

	}
};

struct TestParametersSpecToSmear : public TestParametersGaugeField
{
	int rhoIter;
	double rho;

	TestParametersSpecToSmear(std::vector<double> referenceValueIn, int nsIn, int ntIn, FillType fillTypeIn, int rhoIterIn, double rhoIn) :
		TestParametersGaugeField(referenceValueIn, nsIn, ntIn, fillTypeIn)
	{
		referenceValue = referenceValueIn;
		ns = nsIn;
		nt = ntIn;
		fillType = fillTypeIn;
		rhoIter = rhoIterIn;
		rho = rhoIn;
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
				struct TestParametersGaugeField testParams):
			KernelTester(kernelName, hardwareParameters, kernelParameters, testParams),
			numberOfElements(hardwareParameters.getLatticeVolume() * NDIM)
	{
		gaugefieldBuffer = new hardware::buffers::SU3( numberOfElements, device);
		const Matrixsu3 * gf_host = createGaugefield(numberOfElements, testParams.fillType);
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
							struct TestParametersSpecToPlaq testParams):
				GaugefieldTester("plaquette", hardwareParameters, kernelParameters, testParams), typeOfPlaquette(testParams.typeOfPlaquette)  {
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

	void instantiateMockupsAndCallTester(struct TestParametersSpecToPlaq testParams)
	{
		hardware::HardwareParametersMockup hardwareParameters(testParams.ns,testParams.nt);
		hardware::code::OpenClKernelParametersMockup kernelParameters(testParams.ns,testParams.nt);
		PlaquetteTester plaquetteTester(hardwareParameters, kernelParameters, testParams);
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_1 )
	{
		instantiateMockupsAndCallTester(TestParametersSpecToPlaq {referenceValues{1536.}, ns4, nt4, FillType::cold, 1});
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_2 )
	{
		instantiateMockupsAndCallTester(TestParametersSpecToPlaq {referenceValues{1536.002605}, ns4, nt4, FillType::nonTrivial, 1});
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_TEMPORAL_1 )
	{
		instantiateMockupsAndCallTester(TestParametersSpecToPlaq {referenceValues{768.}, ns4, nt4, FillType::cold, 2});
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_TEMPORAL_2 )
	{
		instantiateMockupsAndCallTester(TestParametersSpecToPlaq {referenceValues{768.00130250240136}, ns4, nt4, FillType::nonTrivial, 2});
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_SPATIAL_1 )
	{
		instantiateMockupsAndCallTester(TestParametersSpecToPlaq {referenceValues{768.}, ns4, nt4, FillType::cold, 3});
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_SPATIAL_2 )
	{
		instantiateMockupsAndCallTester(TestParametersSpecToPlaq {referenceValues{768.00130250240136}, ns4, nt4, FillType::nonTrivial, 3});
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_REDUCTION_1 )
	{
		instantiateMockupsAndCallTester(TestParametersSpecToPlaq {referenceValues{24576}, ns8, nt8, FillType::cold, 1});
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_REDUCTION_2 )
	{
		instantiateMockupsAndCallTester(TestParametersSpecToPlaq {referenceValues{124416}, ns12, nt12, FillType::cold, 1});
	}

	BOOST_AUTO_TEST_CASE( PLAQUETTE_REDUCTION_3 )
	{
		instantiateMockupsAndCallTester(TestParametersSpecToPlaq {referenceValues{393216}, ns16, nt16, FillType::cold, 1});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( POLYAKOV )

	class PolyakovloopTester : public GaugefieldTester {
	public:
			PolyakovloopTester(const hardware::HardwareParametersMockup & hardwareParameters, const hardware::code::OpenClKernelParametersMockup & kernelParameters,
					TestParametersGaugeField testParams):
				GaugefieldTester("polyakov", hardwareParameters, kernelParameters, testParams) {
			const hardware::buffers::Plain<hmc_complex> pol(1, device);
			code->polyakov_device(gaugefieldBuffer, &pol);

			hmc_complex kernelResult_tmp;
			pol.dump(&kernelResult_tmp);
			kernelResult[0] = kernelResult_tmp.re;
			kernelResult[1] = kernelResult_tmp.im;
		}
	};

	void instantiateMockupsAndCallTester(struct TestParametersGaugeField testParams)
	{
		hardware::HardwareParametersMockup hardwareParameters(testParams.ns,testParams.nt);
		hardware::code::OpenClKernelParametersMockup kernelParameters(testParams.ns,testParams.nt);
		PolyakovloopTester polyakovloopTester(hardwareParameters, kernelParameters, testParams);
	}

	BOOST_AUTO_TEST_CASE( POLYAKOV_1 )
	{
		instantiateMockupsAndCallTester(TestParametersGaugeField {referenceValues{64.,0.}, ns4, nt4, FillType::cold});
	}

	BOOST_AUTO_TEST_CASE( POLYAKOV_2)
	{
		instantiateMockupsAndCallTester(TestParametersGaugeField {referenceValues{-17.1117721375,-31.0747993518}, ns4, nt4, FillType::nonTrivial});
	}

	BOOST_AUTO_TEST_CASE( POLYAKOV_REDUCTION_1 )
	{
		instantiateMockupsAndCallTester(TestParametersGaugeField {referenceValues{512.,0.}, ns8, nt8, FillType::cold});
	}

	BOOST_AUTO_TEST_CASE( POLYAKOV_REDUCTION_2 )
	{
		instantiateMockupsAndCallTester(TestParametersGaugeField {referenceValues{1728.,0.}, ns12, nt12, FillType::cold});
	}

	BOOST_AUTO_TEST_CASE( POLYAKOV_REDUCTION_3 )
	{
		instantiateMockupsAndCallTester(TestParametersGaugeField {referenceValues{4096.,0.}, ns16, nt16, FillType::cold});
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
					struct TestParametersSpecToSmear testParams):
				GaugefieldTester("stout_smear", hardwareParameters, kernelParameters, testParams) {

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

void instantiateMockupsAndCallTester(struct TestParametersSpecToSmear testParams)
	{
		const bool useSmearing = true;
		hardware::HardwareParametersMockup hardwareParameters(testParams.ns,testParams.nt);
		hardware::code::OpenClKernelParametersMockup kernelParameters(testParams.ns,testParams.nt, testParams.rhoIter, testParams.rho, useSmearing);
		StoutSmearTester StoutSmearTester(hardwareParameters, kernelParameters, testParams);
	}

	BOOST_AUTO_TEST_CASE( STOUT_SMEAR_1 )
	{
		const int rhoIterIn = 1;
		const double rhoIn = 0.001;
		instantiateMockupsAndCallTester(TestParametersSpecToSmear {referenceValues{-1234}, ns4, nt4, FillType::nonTrivial, rhoIterIn, rhoIn});
	}

	BOOST_AUTO_TEST_CASE( STOUT_SMEAR_2 )
	{
		const int rhoIterIn = 1;
		const double rhoIn = 0.;
		instantiateMockupsAndCallTester(TestParametersSpecToSmear {referenceValues{-1234}, ns4, nt4, FillType::nonTrivial, rhoIterIn, rhoIn});
	}

	BOOST_AUTO_TEST_CASE( STOUT_SMEAR_3 )
	{
		const int rhoIterIn = 1;
		const double rhoIn = 0.001538;
		instantiateMockupsAndCallTester(TestParametersSpecToSmear {referenceValues{-1234}, ns4, nt4, FillType::nonTrivial, rhoIterIn, rhoIn});
	}

	BOOST_AUTO_TEST_CASE( STOUT_SMEAR_4 )
	{
		const int rhoIterIn = 1;
		const double rhoIn = 0.;
		instantiateMockupsAndCallTester(TestParametersSpecToSmear {referenceValues{1536}, ns4, nt4, FillType::cold, rhoIterIn, rhoIn});
	}

	BOOST_AUTO_TEST_CASE( STOUT_SMEAR_5 )
	{
		const int rhoIterIn = 1;
		const double rhoIn = 0.001;
		instantiateMockupsAndCallTester(TestParametersSpecToSmear {referenceValues{1536}, ns4, nt4, FillType::cold, rhoIterIn, rhoIn});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( RECTANGLES )

	class RectanglesTester : public GaugefieldTester {
	public:
			RectanglesTester(const hardware::HardwareParametersMockup & hardwareParameters, const hardware::code::OpenClKernelParametersMockup & kernelParameters,
					struct TestParametersGaugeField testParams):
				GaugefieldTester("rectangles", hardwareParameters, kernelParameters, testParams) {
			const hardware::buffers::Plain<hmc_float> rect(1, device );
			code->rectangles_device(gaugefieldBuffer, &rect);
			rect.dump(&kernelResult[0]);
		}
	};

	void instantiateMockupsAndCallTester(struct TestParametersGaugeField testParams)
	{
		const bool useRectangles = true;
		hardware::HardwareParametersMockup hardwareParameters(testParams.ns,testParams.nt);
		hardware::code::OpenClKernelParametersMockup kernelParameters(testParams.ns,testParams.nt, useRectangles);
		RectanglesTester rectanglesTester(hardwareParameters, kernelParameters, testParams);
	}

	BOOST_AUTO_TEST_CASE( RECTANGLES_1 )
	{
		instantiateMockupsAndCallTester(TestParametersGaugeField {referenceValues{3072.}, ns4, nt4, FillType::cold});
	}

	BOOST_AUTO_TEST_CASE( RECTANGLES_2 )
	{
		instantiateMockupsAndCallTester(TestParametersGaugeField {referenceValues{3072.00781501}, ns4, nt4, FillType::nonTrivial});
	}

	BOOST_AUTO_TEST_CASE( RECTANGLES_REDUCTION_1 )
	{
			instantiateMockupsAndCallTester(TestParametersGaugeField {referenceValues{49152.}, ns8, nt8, FillType::cold});
	}

	BOOST_AUTO_TEST_CASE( RECTANGLES_REDUCTION_2 )
	{
		instantiateMockupsAndCallTester(TestParametersGaugeField {referenceValues{248832.}, ns12, nt12, FillType::cold});
	}

	BOOST_AUTO_TEST_CASE( RECTANGLES_REDUCTION_3 )
	{
		instantiateMockupsAndCallTester(TestParametersGaugeField {referenceValues{786432.}, ns16, nt16, FillType::cold});
	}

BOOST_AUTO_TEST_SUITE_END()

