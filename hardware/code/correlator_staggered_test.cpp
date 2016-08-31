/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra,
 * Francesca Cuteri
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

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE_CORRELATORS_STAGGERED

#include "correlator_staggered.hpp"
#include "SpinorStaggeredTester.hpp"
#include "PrngSpinorTester.hpp"

struct StaggeredFermionsSourceTestParameters : public PrngSpinorStaggeredTestParameters
{
	StaggeredFermionsSourceTestParameters(const LatticeExtents lE, common::sourcecontents sC, const int iterations):
		TestParameters(lE, 10e-2), PrngSpinorStaggeredTestParameters(lE, iterations, false), sourcecontent(sC) {}; // In calling the TestParameters ctor, the testPrecision is reduced, so as related tests can pass with a reasonable number of iterations!
	common::sourcecontents sourcecontent;
};

int countNonZeroElements(const su3vec * in, const int numberOfElements)
{
	int result = 0;
	for (int i = 0; i< numberOfElements; i++)
	{
		if( ! in[i].e0.re == 0.) //su3vec are uniformly filled so the check is done on the first component only
			result += 1;
		if( ! in[i].e0.im == 0.)
			result += 1;
	}
	return result;
}

double normalize(double valueIn, const LatticeExtents lE)
{
	return valueIn/= calculateSpinorfieldSize(lE) * 6;
}

ReferenceValues calculateReferenceValues_volumeSource(const StaggeredFermionsSourceTestParameters & tP)
{
	double mean, variance, nonzeroEntries;
	if (tP.sourcecontent == common::sourcecontents::one)
	{
		mean = normalize(3 * calculateSpinorfieldSize(tP.latticeExtents), tP.latticeExtents);
	}
	else if (tP.sourcecontent == common::sourcecontents::gaussian or tP.sourcecontent == common::sourcecontents::z4 or tP.sourcecontent == common::sourcecontents::z2)
	{
		mean = 0.;
	}
	else
		mean = 0.123456;

	if (tP.sourcecontent == common::sourcecontents::gaussian or tP.sourcecontent == common::sourcecontents::z4)
	{
		nonzeroEntries = calculateEvenOddSpinorfieldSize(tP.latticeExtents) * 2; //real and imaginary parts count separately in the counting
	}
	else if (tP.sourcecontent == common::sourcecontents::z2 or tP.sourcecontent == common::sourcecontents::one)
	{
		nonzeroEntries = calculateEvenOddSpinorfieldSize(tP.latticeExtents);
	}
	else
		nonzeroEntries = 123456;
	variance = sqrt(normalize(((0. - mean) * (0. - mean) * 3 + (1. - mean) * (1. - mean) * 3) * calculateSpinorfieldSize(tP.latticeExtents), tP.latticeExtents));
	return ReferenceValues{mean, variance, nonzeroEntries};
}

struct SourceTester : public PrngSpinorStaggeredTester
{
	SourceTester(const std::string kernelName, const ParameterCollection pC, const StaggeredFermionsSourceTestParameters tP, const ReferenceValues rV):
	     PrngSpinorStaggeredTester(kernelName, pC, tP, calculateEvenOddSpinorfieldSize(tP.latticeExtents), rV), numberOfNonZeroEntries(0)
	{
		code = SpinorStaggeredTester::device->getCorrelatorStaggeredCode();
		outSpinor = new hardware::buffers::SU3vec(tP.latticeExtents, SpinorStaggeredTester::device);
	}
	
	virtual ~SourceTester(){
		numberOfNonZeroEntries = countNonZeroElements (&hostOutput[0], numberOfElements);
		kernelResult.at(2) = numberOfNonZeroEntries;
	}
	
   protected:
	int numberOfNonZeroEntries;
	const hardware::code::Correlator_staggered * code;
	const hardware::buffers::SU3vec *outSpinor;
};

struct VolumeSourceTester : public SourceTester
{
	VolumeSourceTester(const ParameterCollection pC, const StaggeredFermionsSourceTestParameters tP, const int numberOfElements) :
		SourceTester("Volume_source", pC, tP, calculateReferenceValues_volumeSource(tP))
	{
		for (unsigned int i = 0; i< tP.iterations; i++){
		  outSpinor->clear();
		  code->create_volume_source_stagg_eoprec_device(outSpinor, prngStates);
		  outSpinor->dump(&hostOutput[i*numberOfElements]);
		}
	}
};

void testVolumeSource(const LatticeExtents lE, const common::sourcecontents sC, const int iterations)
{
	StaggeredFermionsSourceTestParameters parametersForThisTest(lE, sC, iterations);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.latticeExtents, true);
	hardware::code::OpenClKernelParametersMockupForStaggeredSourceTests kernelParameters(parametersForThisTest.latticeExtents, parametersForThisTest.sourcecontent);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	VolumeSourceTester(parameterCollection, parametersForThisTest, calculateEvenOddSpinorfieldSize(lE));
}


BOOST_AUTO_TEST_SUITE(SRC_VOLUME)

	BOOST_AUTO_TEST_CASE( SRC_VOLUME_1 )
	{
		testVolumeSource(LatticeExtents{ns4, nt4}, common::sourcecontents::one, 500);
	}

	BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(SRC_VOLUME_2, 2)

	BOOST_AUTO_TEST_CASE( SRC_VOLUME_2 )
	{
		testVolumeSource(LatticeExtents{ns4, nt4}, common::sourcecontents::z4, 1000);
	}

	BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(SRC_VOLUME_3, 2)

	BOOST_AUTO_TEST_CASE( SRC_VOLUME_3 )
	{
		testVolumeSource(LatticeExtents{ns4, nt4}, common::sourcecontents::gaussian, 2000);
	}

	BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(SRC_VOLUME_4, 2)

	BOOST_AUTO_TEST_CASE( SRC_VOLUME_4 )
	{
		testVolumeSource(LatticeExtents{ns4, nt4}, common::sourcecontents::z2, 1000);
	}

BOOST_AUTO_TEST_SUITE_END()


