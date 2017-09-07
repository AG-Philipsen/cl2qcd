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

ReferenceValues calculateReferenceValues_squarenorm(const LatticeExtents lE , const GaugeMomentumFilltype fillTypesIn)
{
	switch( fillTypesIn )
	{
		case GaugeMomentumFilltype::One :
		{
			return ReferenceValues{lE.getLatticeVolume() * NDIM * 8.};
		}
		case GaugeMomentumFilltype::Ascending:
		{
			return ReferenceValues{lE.getLatticeVolume() * NDIM * sumOfIntegersSquared(8)};
		}
		default:
		{
			return defaultReferenceValues();
		}
	}
}

ReferenceValues calculateReferenceValues_set_zero()
{
	return ReferenceValues{0.};
}

ReferenceValues calculateReferenceValues_saxpy(const LatticeExtents lE , double alphaIn)
{
	return ReferenceValues{(1. + alphaIn) * (1. + alphaIn) * lE.getLatticeVolume() * NDIM * sumOfIntegersSquared(8)};
}

ReferenceValues calculateReferenceValues_gaussian()
{
	return ReferenceValues{ 0., 1.};
}

struct SquarenormTester : public GaugemomentumTester
{
  SquarenormTester(const ParameterCollection pC, const GaugemomentumTestParameters tP) :
  GaugemomentumTester("gaugemomenta squarenorm", pC, calculateReferenceValues_squarenorm(tP.latticeExtents, tP.fillType), tP)
  {
	GaugemomentumCreator gm(tP.latticeExtents);
    gaugemomentumBuffer = new hardware::buffers::Gaugemomentum(tP.latticeExtents, this->device);
    hardware::buffers::Gaugemomentum in(tP.latticeExtents, device);
    code->importGaugemomentumBuffer(gaugemomentumBuffer, reinterpret_cast<ae*>( gm.createGaugemomentumBasedOnFilltype(tP.fillType) ));
    calcSquarenormAndStoreAsKernelResult(gaugemomentumBuffer);
  }
private:
    hardware::buffers::Gaugemomentum * gaugemomentumBuffer;
};

struct SetZeroTester : public GaugemomentumTester
{
  SetZeroTester(const ParameterCollection pC, const GaugemomentumTestParameters tP) :
  GaugemomentumTester("set zero", pC, calculateReferenceValues_set_zero(), tP)
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
  SaxpyTester(const ParameterCollection pC, const GaugemomentumTestParameters tP) :
  GaugemomentumTester("saxpy", pC, calculateReferenceValues_saxpy(tP.latticeExtents, tP.coefficient), tP)
  {
	GaugemomentumCreator gm(tP.latticeExtents);
	gaugemomentumBuffer = new hardware::buffers::Gaugemomentum(tP.latticeExtents, this->device);
    hardware::buffers::Gaugemomentum out(tP.latticeExtents, device);
    code->importGaugemomentumBuffer(gaugemomentumBuffer, reinterpret_cast<ae*>( gm.createGaugemomentumBasedOnFilltype(tP.fillType) ));
    code->importGaugemomentumBuffer(&out, reinterpret_cast<ae*>( gm.createGaugemomentumBasedOnFilltype(tP.fillType) ));
    doubleBuffer->load(&tP.coefficient);

    code->saxpy_device(gaugemomentumBuffer, &out, doubleBuffer, &out);
    calcSquarenormAndStoreAsKernelResult(&out);
  }
private:
  hardware::buffers::Gaugemomentum * gaugemomentumBuffer;
};

struct PrngGaugemomentumTestParameters: public GaugemomentumTestParameters
{
	PrngGaugemomentumTestParameters(const LatticeExtents lE, const int iterationsIn) :
		GaugemomentumTestParameters(lE, 10e-4), iterations(iterationsIn) {};
	const unsigned int iterations;
};

struct PrngGaugemomentumTester : public GaugemomentumTester
{
	PrngGaugemomentumTester(const std::string kernelName, const ParameterCollection pC, const PrngGaugemomentumTestParameters &tP, const int numberOfElements):
		GaugemomentumTester(kernelName, pC, calculateReferenceValues_gaussian(), tP),
		numberOfElements(numberOfElements),
		mean(0.), variance(0.),
		hostOutput(std::vector<ae> (numberOfElements * tP.iterations)),
		testParameters(tP),
		hostSeed(pC.kernelParameters.getHostSeed()),
		useSameRandomNumbers(pC.hardwareParameters.useSameRandomNumbers())

		{
		prng_init(hostSeed);
		prngStates = new hardware::buffers::PRNGBuffer(device, useSameRandomNumbers );
		auto codePrng = device->getPrngCode();
		codePrng->initialize(prngStates, hostSeed);
		}

		~PrngGaugemomentumTester()
		{
			calculateMean();
			calculateVariance();

			kernelResult.at(0) = mean;
			kernelResult.at(1) = sqrt(variance);

			delete prngStates;
		}

		double normalize (double valueIn)
		{
			return valueIn /= testParameters.iterations * numberOfElements * 8.;
		}

		void calculateMean()
		{
			for (unsigned int i = 0; i < testParameters.iterations; i++) {
				mean += count_gm(&hostOutput[i * numberOfElements], numberOfElements);
			}
			mean = normalize(mean);
		}

		void calculateVariance()
		{
			for (unsigned int i = 0; i < testParameters.iterations; i++) {
				variance += calc_var_gm(&hostOutput[i * numberOfElements], numberOfElements, mean);
			}
			variance = normalize(variance);
		}
protected:
	const int numberOfElements;
	double mean, variance;
	std::vector<ae> hostOutput;
	const PrngGaugemomentumTestParameters & testParameters;
	const hardware::buffers::PRNGBuffer* prngStates;
private:
	uint32_t hostSeed;
	bool useSameRandomNumbers;
};

struct GaussianTester : public PrngGaugemomentumTester
{
  GaussianTester(const ParameterCollection pC, const PrngGaugemomentumTestParameters tP) :
  PrngGaugemomentumTester("gaussian gaugemomentum", pC, tP, calculateGaugemomentumSize(tP.latticeExtents)){}
  ~GaussianTester()
  {
	const hardware::buffers::Gaugemomentum gm_out(numberOfElements, device);
	for (unsigned int i = 0; i< testParameters.iterations; i++){
	  code->generate_gaussian_gaugemomenta_device(&gm_out, prngStates);
	  gm_out.dump(&hostOutput[i*numberOfElements]);
	}
  }
};

template<class TesterClass>
void callTest(const LatticeExtents lE)
{
	GaugemomentumTestParameters parametersForThisTest(lE);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.latticeExtents);
	hardware::code::OpenClKernelParametersMockupForSpinorTests kernelParameters(parametersForThisTest.latticeExtents);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, parametersForThisTest);
}

template<class TesterClass>
void performTest(const LatticeExtents lE, const int iterations)
{
	PrngGaugemomentumTestParameters parametersForThisTest(lE, iterations);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.latticeExtents);
	hardware::code::OpenClKernelParametersMockupForSpinorTests kernelParameters(parametersForThisTest.latticeExtents);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, parametersForThisTest);
}

template<class TesterClass>
void callTest(const LatticeExtents lE, const GaugeMomentumFilltype gmF)
{
	GaugemomentumTestParameters parametersForThisTest(lE, gmF);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.latticeExtents);
	hardware::code::OpenClKernelParametersMockupForSpinorTests kernelParameters(parametersForThisTest.latticeExtents);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, parametersForThisTest);
}

template<class TesterClass>
void callTest(const LatticeExtents lE, const GaugeMomentumFilltype gmF, const double alphaIn)
{
	GaugemomentumTestParameters parametersForThisTest(lE, gmF, alphaIn);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.latticeExtents);
	hardware::code::OpenClKernelParametersMockupForSpinorTests kernelParameters(parametersForThisTest.latticeExtents);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, parametersForThisTest);
}

void testSquarenorm(const LatticeExtents lE, const GaugeMomentumFilltype gmF)
{
	callTest<SquarenormTester>(lE, gmF);
}

void testSetZero(const LatticeExtents lE)
{
	callTest<SetZeroTester>(lE);
}

void testSaxpy(const LatticeExtents lE, const GaugeMomentumFilltype gmF, const double alpha)
{
	callTest<SaxpyTester>(lE, gmF, alpha);
}

void testGaussianGaugemomentum(const LatticeExtents lE, const int iterations)
{
	performTest<GaussianTester>(lE, iterations);
}

BOOST_AUTO_TEST_SUITE( SQUARENORM )

	BOOST_AUTO_TEST_CASE(SQUARENORM_1  )
	{
		testSquarenorm(LatticeExtents{ns4, nt4}, GaugeMomentumFilltype::One);
	}

	BOOST_AUTO_TEST_CASE(SQUARENORM_2  )
	{
		testSquarenorm(LatticeExtents{ns4, nt4}, GaugeMomentumFilltype::Ascending);
	}
	BOOST_AUTO_TEST_CASE(SQUARENORM_REDUCTION_1  )
	{
		testSquarenorm(LatticeExtents{ns16, nt8}, GaugeMomentumFilltype::One);
	}

	BOOST_AUTO_TEST_CASE(SQUARENORM_REDUCTION_2  )
	{
		testSquarenorm(LatticeExtents{ns8, nt4}, GaugeMomentumFilltype::Ascending);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( SET_ZERO )

	BOOST_AUTO_TEST_CASE( SET_ZERO_1 )
	{
		testSetZero(LatticeExtents{ns4, nt4});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( SAXPY )

	BOOST_AUTO_TEST_CASE( SAXPY_1 )
	{
		testSaxpy(LatticeExtents{ns4, nt8}, GaugeMomentumFilltype::Ascending, nonTrivialParameter);
	}

	BOOST_AUTO_TEST_CASE( SAXPY_2 )
	{
		testSaxpy(LatticeExtents{ns8, nt4}, GaugeMomentumFilltype::Ascending, -nonTrivialParameter);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(GENERATE_GAUSSIAN_GAUGEMOMENTA  )

BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(GENERATE_GAUSSIAN_GAUGEMOMENTA_1, 2)

	BOOST_AUTO_TEST_CASE(GENERATE_GAUSSIAN_GAUGEMOMENTA_1 )
	{
	  testGaussianGaugemomentum(LatticeExtents{ns8, nt4}, 1000);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( CONVERT_TO_SOA )

	BOOST_AUTO_TEST_CASE( CONVERT_TO_SOA_1 )
	{
		BOOST_TEST_MESSAGE("NOT YET IMPLEMENTED!!");
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( CONVERT_FROM_SOA )

	BOOST_AUTO_TEST_CASE( CONVERT_FROM_SOA_1 )
	{
		BOOST_TEST_MESSAGE("NOT YET IMPLEMENTED!!");
	}

BOOST_AUTO_TEST_SUITE_END()

