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

enum CorrelatorDirection {temporal=0, spatialX=1, spatialY=2, spatialZ=3};
typedef std::string KernelIdentifier;


struct StaggeredFermionsSourceTestParameters : public PrngSpinorStaggeredTestParameters
{
	StaggeredFermionsSourceTestParameters(const LatticeExtents lE, common::sourcecontents sC, common::sourcetypes sT, const int iterations):
		TestParameters(lE, 10e-2), PrngSpinorStaggeredTestParameters(lE, iterations, false), sourcecontent(sC), sourcetype(sT) {}; // In calling the TestParameters ctor, the testPrecision is reduced, so as related tests can pass with a reasonable number of iterations!
	const common::sourcecontents sourcecontent;
	const common::sourcetypes sourcetype;
};

struct CorrelatorStaggeredTestParameters : public SpinorStaggeredTestParameters
{
	CorrelatorStaggeredTestParameters(LatticeExtents lE, CorrelatorDirection directionIn, SpinorFillTypes sF) :
		TestParameters(lE), SpinorStaggeredTestParameters(lE, sF), direction(directionIn) {};
	const double kappa=1.;
	CorrelatorDirection direction;
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

double normalizeEvenOdd(double valueIn, const LatticeExtents lE)
{
    return valueIn/= calculateEvenOddSpinorfieldSize(lE) * 6;
}

ReferenceValues calculateReferenceValues_volumeSource(const StaggeredFermionsSourceTestParameters tP)
{
	double mean, variance, nonzeroEntries;
	if (tP.sourcecontent == common::sourcecontents::one)
	{
		mean = normalizeEvenOdd(3 * calculateSpinorfieldSize(tP.latticeExtents), tP.latticeExtents);
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
	variance = sqrt(normalizeEvenOdd(((0. - mean) * (0. - mean) * 3 + (1. - mean) * (1. - mean) * 3) * calculateSpinorfieldSize(tP.latticeExtents), tP.latticeExtents));
	return ReferenceValues{mean, variance, nonzeroEntries};
}

ReferenceValues calculateReferenceValues_pointSource(const StaggeredFermionsSourceTestParameters tP)
{
	double mean, variance;
	if (tP.sourcecontent == common::sourcecontents::one)
	{
		mean = normalizeEvenOdd(1, tP.latticeExtents);
	}
	else
		mean = 0.123456;

	variance = sqrt((0. - mean) * (0. - mean) * ((calculateEvenOddSpinorfieldSize(tP.latticeExtents) - 1) * 24 + 23));

	return ReferenceValues{mean, variance, 1};
}


struct SourceTester : public PrngSpinorStaggeredTester
{
	SourceTester(const std::string kernelName, const ParameterCollection pC, const StaggeredFermionsSourceTestParameters & tP, const ReferenceValues rV):
	     PrngSpinorStaggeredTester(kernelName, pC, tP, calculateEvenOddSpinorfieldSize(tP.latticeExtents), rV), numberOfNonZeroEntries(0)
	{
		code = SpinorStaggeredTester::device->getCorrelatorStaggeredCode();
		outSpinor = new hardware::buffers::SU3vec(tP.latticeExtents, SpinorStaggeredTester::device);
	}
	
	virtual ~SourceTester(){
		numberOfNonZeroEntries = countNonZeroElements(&hostOutput[0], numberOfElements);
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

struct PointSourceTester : public SourceTester
{
	PointSourceTester(const ParameterCollection pC, const StaggeredFermionsSourceTestParameters tP, const int numberOfElements):
		SourceTester("Point_source", pC, tP, calculateReferenceValues_pointSource(tP))
	{
		for (unsigned int i = 0; i < testParameters.iterations; i++)
		{
			code->create_point_source_stagg_eoprec_device(outSpinor, i, 0, 0);
			outSpinor->dump(&hostOutput[i * numberOfElements]);
		}
	}
};

void testPointSource(const LatticeExtents lE, const common::sourcecontents sC,  const int iterations)
{
	StaggeredFermionsSourceTestParameters parametersForThisTest(lE, sC, common::sourcetypes::point, iterations);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.latticeExtents, true);
	hardware::code::OpenClKernelParametersMockupForStaggeredSourceTests kernelParameters(parametersForThisTest.latticeExtents, parametersForThisTest.sourcecontent, common::sourcetypes::point);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	PointSourceTester(parameterCollection, parametersForThisTest, calculateEvenOddSpinorfieldSize(lE));
}


void testVolumeSource(const LatticeExtents lE, const common::sourcecontents sC, const int iterations)
{
	StaggeredFermionsSourceTestParameters parametersForThisTest(lE, sC, common::sourcetypes::volume, iterations);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.latticeExtents, true);
	hardware::code::OpenClKernelParametersMockupForStaggeredSourceTests kernelParameters(parametersForThisTest.latticeExtents, parametersForThisTest.sourcecontent, common::sourcetypes::volume);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	VolumeSourceTester(parameterCollection, parametersForThisTest, calculateEvenOddSpinorfieldSize(lE));
}


struct CorrelatorStaggeredTester : public EvenOddSpinorStaggeredTester
{
	CorrelatorStaggeredTester(const KernelIdentifier kI, const ParameterCollection pC, const ReferenceValues rV, const CorrelatorStaggeredTestParameters tP):
		EvenOddSpinorStaggeredTester(kI, pC, tP, rV), correlatorEntries(rV.size())
	{
		code = device->getCorrelatorStaggeredCode();
		CorrelatorResult = new hardware::buffers::Plain<hmc_float>(correlatorEntries, device);
		CorrelatorResult->clear();

		for( unsigned int i=0; i<tP.fillTypes.size(); i++ )
		{
			EvenOddSpinorStaggeredfieldCreator sf(tP.latticeExtents);
			spinorStaggeredfields.push_back( new hardware::buffers::SU3vec (sf.numberOfElements, device) );
			auto spinorStaggeredfield = sf.createSpinorfield(tP.fillTypes.at(i));
			spinorStaggeredfields.at(i)->load(spinorStaggeredfield);
			delete[] spinorStaggeredfield;
		}

		code->pseudoScalarCorrelator(CorrelatorResult, spinorStaggeredfields.at(0));

	};
	~CorrelatorStaggeredTester()
	{
		CorrelatorResult->dump(&kernelResult.at(0));
	}
protected:
	const int correlatorEntries;
	const hardware::buffers::Plain<hmc_float> * CorrelatorResult;
	const hardware::code::Correlator_staggered * code;
	std::vector< const hardware::buffers::SU3vec* > spinorStaggeredfields;
};

//struct ComponentwiseCorrelatorStaggeredTester : public CorrelatorStaggeredTester
//{
//	ComponentwiseCorrelatorStaggeredTester(const KernelIdentifier kI, const ParameterCollection pC, const ReferenceValues rV, const CorrelatorStaggeredTestParameters tP):
//		CorrelatorStaggeredTester(kI, pC, rV, tP)
//	{
//		code->pseudoScalarCorrelator(CorrelatorResult, spinorStaggeredfields.at(0) );
//	}
//};


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


BOOST_AUTO_TEST_SUITE(SRC_POINT)

	BOOST_AUTO_TEST_CASE( SRC_POINT_1 )
	{
		testPointSource( LatticeExtents{ns4, nt4}, common::sourcecontents::one, 3);
	}

BOOST_AUTO_TEST_SUITE_END()

