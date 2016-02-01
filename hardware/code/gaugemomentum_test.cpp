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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE HARDWARE_CODE_GAUGEMOMENTUM

#include "GaugemomentumTester.hpp"

struct SquarenormTester : public GaugemomentumTester
{
  SquarenormTester(const ParameterCollection pC, const ReferenceValues rV, const GaugemomentumTestParameters tP) :
    GaugemomentumTester("gaugemomenta squarenorm", pC, rV, tP)
  {
	GaugemomentumCreator gm(tP.latticeExtents);
    gaugemomentumBuffer = new hardware::buffers::Gaugemomentum(tP.latticeExtents, this->device);
    hardware::buffers::Gaugemomentum in(tP.latticeExtents, device);
    code->importGaugemomentumBuffer(gaugemomentumBuffer, reinterpret_cast<ae*>( gm.createGaugemomentumBasedOnFilltype(GaugeMomentumFilltype::One) ));
    calcSquarenormAndStoreAsKernelResult(gaugemomentumBuffer);
  }
private:
    hardware::buffers::Gaugemomentum * gaugemomentumBuffer;
};

struct SetZeroTester : public GaugemomentumTester
{
  SetZeroTester(const ParameterCollection pC, const ReferenceValues rV, const GaugemomentumTestParameters tP) :
    GaugemomentumTester("set zero", pC, rV, tP)
  {
	GaugemomentumCreator gm(tP.latticeExtents);
	gaugemomentumBuffer = new hardware::buffers::Gaugemomentum(tP.latticeExtents, this->device);
    code->importGaugemomentumBuffer(gaugemomentumBuffer, reinterpret_cast<ae*>( gm.createGaugemomentumBasedOnFilltype(GaugeMomentumFilltype::One) ));
    code->set_zero_gaugemomentum(gaugemomentumBuffer);
    calcSquarenormAndStoreAsKernelResult(gaugemomentumBuffer);
  }
private:
  hardware::buffers::Gaugemomentum * gaugemomentumBuffer;
};

struct SaxpyTester : public GaugemomentumTester
{
  SaxpyTester(const ParameterCollection pC, const ReferenceValues rV, const GaugemomentumTestParameters tP) :
    GaugemomentumTester("saxpy", pC, rV, tP)
  {
	GaugemomentumCreator gm(tP.latticeExtents);
	gaugemomentumBuffer = new hardware::buffers::Gaugemomentum(tP.latticeExtents, this->device);
    hardware::buffers::Gaugemomentum out(tP.latticeExtents, device);
    code->importGaugemomentumBuffer(gaugemomentumBuffer, reinterpret_cast<ae*>( gm.createGaugemomentumBasedOnFilltype(GaugeMomentumFilltype::One) ));
    code->importGaugemomentumBuffer(&out, reinterpret_cast<ae*>( gm.createGaugemomentumBasedOnFilltype(GaugeMomentumFilltype::One) ));
    double alpha = 1.2345; //@todo: make adjustable
    doubleBuffer->load(&alpha);

    code->saxpy_device(gaugemomentumBuffer, &out, doubleBuffer, &out);
    calcSquarenormAndStoreAsKernelResult(&out);
  }
private:
  hardware::buffers::Gaugemomentum * gaugemomentumBuffer;
};

struct GaussianTester : public GaugemomentumTester
{
  GaussianTester(const ParameterCollection pC, const ReferenceValues rV, const GaugemomentumTestParameters tP) :
    GaugemomentumTester("gaussian gaugemomentum", pC, rV, tP), hostSeed(1234), useSameRandomNumbers(false)
  {
	//@todo: This must more similar to PrngSpinorTester! I only did quick modifications.

	GaugemomentumCreator gm(tP.latticeExtents);
	gaugemomentumBuffer = new hardware::buffers::Gaugemomentum(tP.latticeExtents, this->device);
	prng_init(hostSeed);
	prngStates = new hardware::buffers::PRNGBuffer(device, useSameRandomNumbers );
	auto codePrng = device->getPrngCode();
	codePrng->initialize(prngStates, hostSeed);

	double result = 0.;
	double sum = 0.;
	int iterations = 10; //todo: make more general
	ae * gm_out;

	int numberOfGaugemomentumElements = calculateGaugemomentumSize(tP.latticeExtents);
	gm_out = new ae[numberOfGaugemomentumElements * iterations];
	BOOST_REQUIRE(gm_out);

	for (int i = 0; i< iterations; i++){
	  code->generate_gaussian_gaugemomenta_device(gaugemomentumBuffer, prngStates);
	  gaugemomentumBuffer->dump(&gm_out[i*numberOfGaugemomentumElements]);
	  sum += gm.count_gm(&gm_out[i*numberOfGaugemomentumElements], numberOfGaugemomentumElements);
	}
	sum = sum/iterations/numberOfGaugemomentumElements/8;
	result= sum;

	  double var=0.;
	  for (int i=0; i<iterations; i++){
	    var += gm.calc_var_gm(&gm_out[i*numberOfGaugemomentumElements], numberOfGaugemomentumElements, sum);
	  }
	  var=var/iterations/numberOfGaugemomentumElements/8;

	  result = sqrt(var);

	kernelResult[0] = result;
  }
private:
	const hardware::buffers::PRNGBuffer* prngStates;
	uint32_t hostSeed;
	bool useSameRandomNumbers;
	hardware::buffers::Gaugemomentum * gaugemomentumBuffer;
};

template<class TesterClass>
void callTest(const LatticeExtents lE)
{
	GaugemomentumTestParameters parametersForThisTest(lE);
	//todo: Work over these!
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt, true);
	hardware::code::OpenClKernelParametersMockupForSpinorTests kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt, true);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, defaultReferenceValues(), parametersForThisTest);
}

void testSquarenorm(const LatticeExtents lE)
{
	callTest<SquarenormTester>(lE);
}

void testSetZero(const LatticeExtents lE)
{
	callTest<SetZeroTester>(lE);
}

void testSaxpy(const LatticeExtents lE)
{
	callTest<SaxpyTester>(lE);
}

void testGaussianGaugemomentum(const LatticeExtents lE)
{
	callTest<GaussianTester>(lE);
}

BOOST_AUTO_TEST_SUITE( SQUARENORM )

	//@todo: add non-trivial tests like in "squarenorm_input_1", "squarenorm_input_2", "squarenorm_reduction_input_1", "squarenorm_reduction_input_2", "squarenorm_reduction_input_3"
	BOOST_AUTO_TEST_CASE(SQUARENORM_1  )
	{
		testSquarenorm(LatticeExtents{ns4, nt4});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( SET_ZERO )

	BOOST_AUTO_TEST_CASE( SET_ZERO_1 )
	{
		testSetZero(LatticeExtents{ns4, nt4});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( SAXPY )

	//@todo: add tests like "saxpy_input_{1-5}"
	BOOST_AUTO_TEST_CASE( SAXPY_1 )
	{
		testSaxpy(LatticeExtents{ns4, nt8});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(GENERATE_GAUSSIAN_GAUGEMOMENTA  )

	//@todo: add tests like in "gaussian_input_{1-4}"
	BOOST_AUTO_TEST_CASE(GENERATE_GAUSSIAN_GAUGEMOMENTA_1 )
	{
	  testGaussianGaugemomentum(LatticeExtents{ns8, nt4});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( CONVERT_TO_SOA )

	BOOST_AUTO_TEST_CASE( CONVERT_TO_SOA_1 )
	{
		BOOST_MESSAGE("NOT YET IMPLEMENTED!!");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( CONVERT_FROM_SOA )

	BOOST_AUTO_TEST_CASE( CONVERT_FROM_SOA_1 )
	{
		BOOST_MESSAGE("NOT YET IMPLEMENTED!!");
	}

BOOST_AUTO_TEST_SUITE_END()

