/*
 * Copyright 2012, 2013, 2015 Lars Zeidlewicz, Christopher Pinke,
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

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE_CORRELATORS
#include <boost/test/unit_test.hpp>

#include "correlator.hpp"
#include "PrngSpinorTester.hpp"

enum CorrelatorDirection {temporal=0, spatialX=1, spatialY=2, spatialZ=3}; //@todo: should this be name more general, "direction" or so, and be moved to common?
typedef std::string KernelIdentifier;

struct SourceTestParameters : public PrngSpinorTestParameters
{
	SourceTestParameters(const LatticeExtents lE, common::sourcecontents sC, common::sourcetypes sT, const int iterations) :
		TestParameters(lE, 10e-2), PrngSpinorTestParameters(lE, iterations), sC(sC), sT(sT) {}
	const common::sourcecontents sC;
	const common::sourcetypes sT;
};

struct CorrelatorTestParameters : public SpinorTestParameters
{
	CorrelatorTestParameters(LatticeExtents lE, CorrelatorDirection directionIn, SpinorFillTypes sF) :
		TestParameters(lE), SpinorTestParameters(lE, sF), direction(directionIn) {};
	const double kappa=1.;
	CorrelatorDirection direction;
};

bool compareToZero_su3vec(const su3vec in)
{
    if (in.e0.re == 0. && in.e0.im == 0.)
        if(in.e1.re == 0. && in.e1.im == 0.)
            if(in.e2.re == 0. && in.e2.im == 0.)
                return true;
    return false;
}

bool compareToZero(const spinor in)
{
	if(compareToZero_su3vec(in.e0))
		if(compareToZero_su3vec(in.e1))
			if(compareToZero_su3vec(in.e2))
				if(compareToZero_su3vec(in.e3))
					return true;
	return false;
}

double normalize(double valueIn, const LatticeExtents lE)
{
	return valueIn/= calculateSpinorfieldSize(lE) * 24;
}

int countNonZeroElements(const spinor * in, const int numberOfElements)
{
	int result = 0;
	for (int i = 0; i< numberOfElements; i++)
	{
		if( !compareToZero(in[i]) )
			result += 1;
	}
	return result;
}

ReferenceValues calculateReferenceValues_volumeSource(const SourceTestParameters & tP)
{
	double mean, variance;
	if (tP.sC == common::sourcecontents::one)
	{
		mean = normalize(12 * calculateSpinorfieldSize(tP.latticeExtents), tP.latticeExtents);
	}
	else if (tP.sC == common::sourcecontents::gaussian or tP.sC == common::sourcecontents::z4)
	{
		mean = 0.;
	}
	else
		mean = 0.123456;
	variance = sqrt(normalize(((0. - mean) * (0. - mean) * 12 + (1. - mean) * (1. - mean) * 12) * calculateSpinorfieldSize(tP.latticeExtents), tP.latticeExtents));
	return ReferenceValues{mean, variance, calculateSpinorfieldSize(tP.latticeExtents)};
}

ReferenceValues calculateReferenceValues_zSliceSource(const SourceTestParameters & tP)
{
	double mean, variance;
	if (tP.sC == common::sourcecontents::one)
	{
		mean = normalize(12*tP.latticeExtents.getNs()*tP.latticeExtents.getNs()*tP.latticeExtents.getNt(), tP.latticeExtents);
	}
	else if (tP.sC == common::sourcecontents::gaussian or tP.sC == common::sourcecontents::z4)
	{
		mean = 0.;
	}
	else
		mean = 0.123456;
	variance = sqrt(normalize((0. - mean) * (0. - mean) * (tP.latticeExtents.getNs()*tP.latticeExtents.getNs()*tP.latticeExtents.getNt() * 12)
			+ (1. - mean) * (1. - mean) * (tP.latticeExtents.getNs()*tP.latticeExtents.getNs()*tP.latticeExtents.getNt() * 12)
			+ (0. - mean) * (0. - mean) * (calculateSpinorfieldSize(tP.latticeExtents) - tP.latticeExtents.getNs()*tP.latticeExtents.getNs()*tP.latticeExtents.getNt()) * 24, tP.latticeExtents));
	return ReferenceValues{mean, variance, (double) tP.latticeExtents.getNs()*tP.latticeExtents.getNs()*tP.latticeExtents.getNt()};
}

ReferenceValues calculateReferenceValues_timeSliceSource(const SourceTestParameters & tP)
{
	double mean, variance;
	if (tP.sC == common::sourcecontents::one)
	{
		mean = normalize(12*tP.latticeExtents.getNs()*tP.latticeExtents.getNs()*tP.latticeExtents.getNs(), tP.latticeExtents);
	}
	else if (tP.sC == common::sourcecontents::gaussian or tP.sC == common::sourcecontents::z4)
	{
		mean = 0.;
	}
	else
		mean = 0.123456;
	variance = sqrt(normalize((0. - mean) * (0. - mean) * (tP.latticeExtents.getNs()*tP.latticeExtents.getNs()*tP.latticeExtents.getNs() * 12)
			+ (1. - mean) * (1. - mean) * (tP.latticeExtents.getNs()*tP.latticeExtents.getNs()*tP.latticeExtents.getNs() * 12)
			+ (0. - mean) * (0. - mean) * (calculateSpinorfieldSize(tP.latticeExtents) - tP.latticeExtents.getNs()*tP.latticeExtents.getNs()*tP.latticeExtents.getNt()) * 24, tP.latticeExtents));
	return ReferenceValues{mean, variance, (double) tP.latticeExtents.getNs()*tP.latticeExtents.getNs()*tP.latticeExtents.getNs()};
}

ReferenceValues calculateReferenceValues_pointSource(const SourceTestParameters & tP)
{
	double mean, variance;
	if (tP.sC == common::sourcecontents::one)
	{
		mean = normalize(1, tP.latticeExtents);
		variance = sqrt((0. - mean) * (0. - mean) * ((calculateSpinorfieldSize(tP.latticeExtents) - 1) * 24 + 23));
	}
	else
		mean = 0.123456;
	return ReferenceValues{mean, variance, 1};
}

struct SourceTester : public PrngSpinorTester
{
	SourceTester(KernelIdentifier kI, const ParameterCollection pC, const SourceTestParameters & tP, const int numberOfElements, const ReferenceValues rV):
		PrngSpinorTester(kI, pC, tP, calculateSpinorfieldSize(tP.latticeExtents), rV ), numberOfNonZeroEntries(0)
	{
		code = device->getCorrelatorCode();
		outSpinor = new hardware::buffers::Plain<spinor>(numberOfElements, device);
	}
	~SourceTester()
	{
		numberOfNonZeroEntries = countNonZeroElements(&hostOutput[0], testParameters.iterations*numberOfElements)/testParameters.iterations;
		kernelResult.at(2) = numberOfNonZeroEntries;
	}
protected:
	int numberOfNonZeroEntries;
	const hardware::code::Correlator * code;
	const hardware::buffers::Plain<spinor> * outSpinor;
};

struct VolumeSourceTester : public SourceTester
{
	VolumeSourceTester(const ParameterCollection pC, const SourceTestParameters & tP, const int numberOfElements):
		SourceTester("create_volume_source", pC, tP, calculateSpinorfieldSize(tP.latticeExtents), calculateReferenceValues_volumeSource(tP))
	{
		for (unsigned int i = 0; i < testParameters.iterations; i++)
		{
			code->create_volume_source_device(outSpinor, prngStates);
			outSpinor->dump(&hostOutput[i * numberOfElements]);
		}
	}
};

struct ZSliceSourceTester : public SourceTester
{
	ZSliceSourceTester(const ParameterCollection pC, const SourceTestParameters & tP, const int numberOfElements):
		SourceTester("create_zslice_source", pC, tP, calculateSpinorfieldSize(tP.latticeExtents), calculateReferenceValues_zSliceSource(tP))
	{
		for (unsigned int i = 0; i < testParameters.iterations; i++)
		{
			code->create_zslice_source_device(outSpinor, prngStates, 0); //@todo: this must be chooseable!
			outSpinor->dump(&hostOutput[i * numberOfElements]);
		}
	}
};

struct TimeSliceSourceTester : public SourceTester
{
	TimeSliceSourceTester(const ParameterCollection pC, const SourceTestParameters & tP, const int numberOfElements):
		SourceTester("create_timeslice_source", pC, tP, calculateSpinorfieldSize(tP.latticeExtents), calculateReferenceValues_timeSliceSource(tP))
	{
		for (unsigned int i = 0; i < testParameters.iterations; i++)
		{
			code->create_timeslice_source_device(outSpinor, prngStates, 0); //@todo: this must be chooseable!
			outSpinor->dump(&hostOutput[i * numberOfElements]);
		}
	}
};

struct PointSourceTester : public SourceTester
{
	PointSourceTester(const ParameterCollection pC, const SourceTestParameters & tP, const int numberOfElements):
		SourceTester("create_point_source", pC, tP, calculateSpinorfieldSize(tP.latticeExtents), calculateReferenceValues_pointSource(tP))
	{
		for (unsigned int i = 0; i < testParameters.iterations; i++)
		{
			code->create_point_source_device(outSpinor, i,0,0); //@todo: this must be chooseable!
			outSpinor->dump(&hostOutput[i * numberOfElements]);
		}
	}
};

struct CorrelatorTester : public NonEvenOddSpinorTester
{
	CorrelatorTester(const KernelIdentifier kI, const ParameterCollection pC, const ReferenceValues rV, const CorrelatorTestParameters tP):
		NonEvenOddSpinorTester(kI, pC, tP, rV), correlatorEntries(rV.size())
	{
		code = device->getCorrelatorCode();
		result = new hardware::buffers::Plain<hmc_float>(correlatorEntries, device);
		result->clear();

		for( unsigned int i=0; i<tP.fillTypes.size(); i++ )
		{
			NonEvenOddSpinorfieldCreator sf(tP.latticeExtents);
			spinorfields.push_back( new hardware::buffers::Plain<spinor> (sf.numberOfElements, device) );
			auto spinorfield = sf.createSpinorfield(tP.fillTypes.at(i));
			spinorfields.at(i)->load(spinorfield);
			delete[] spinorfield;
		}
	};
	~CorrelatorTester()
	{
		result->dump(&kernelResult.at(0));
	}
protected:
	const int correlatorEntries;
	const hardware::buffers::Plain<hmc_float> * result;
	const hardware::code::Correlator * code;
	std::vector< const hardware::buffers::Plain<spinor>* > spinorfields;
};

struct ComponentwiseCorrelatorTester : public CorrelatorTester
{
	ComponentwiseCorrelatorTester(const KernelIdentifier kI, const ParameterCollection pC, const ReferenceValues rV, const CorrelatorTestParameters tP):
		CorrelatorTester(kI, pC, rV, tP)
	{
		code->correlator(code->get_correlator_kernel(kI), result, spinorfields.at(0) );
	}
};

struct ColorwiseCorrelatorTester : public CorrelatorTester
{
	ColorwiseCorrelatorTester(const KernelIdentifier kI, const ParameterCollection pC, const ReferenceValues rV, const CorrelatorTestParameters tP):
		CorrelatorTester(kI, pC, rV, tP)
	{
		code->correlator(code->get_correlator_kernel(kI), result, spinorfields.at(0), spinorfields.at(1), spinorfields.at(2), spinorfields.at(3));
	}
};

template <class TesterClass>
void callTest(const KernelIdentifier kI, const LatticeExtents lE, const CorrelatorDirection cD, const SpinorFillTypes sF, const ReferenceValues rV)
{
	CorrelatorTestParameters parametersForThisTest(lE, cD, sF);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt, false);
	hardware::code::OpenClKernelParametersMockupForCorrelators kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt, parametersForThisTest.kappa, parametersForThisTest.direction);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(kI, parameterCollection, rV, parametersForThisTest );
}

template < class TesterClass >
void callTest(const LatticeExtents lE, const common::sourcecontents sC, const common::sourcetypes sT, const int iterations)
{
	SourceTestParameters parametersForThisTest(lE, sC, sT, iterations);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt, false);
	hardware::code::OpenClKernelParametersMockupForSourceTests kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt, parametersForThisTest.sC, parametersForThisTest.sT);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, parametersForThisTest, calculateSpinorfieldSize(lE));
}

void testVolumeSource(const LatticeExtents lE, const common::sourcecontents sC, const int iterations)
{
	callTest<VolumeSourceTester>( lE, sC, common::sourcetypes::volume, iterations);
}

void testZSliceSource(const LatticeExtents lE, const common::sourcecontents sC, const int iterations)
{
	callTest<ZSliceSourceTester>( lE, sC, common::sourcetypes::zslice, iterations);
}

void testTimeSliceSource(const LatticeExtents lE, const common::sourcecontents sC, const int iterations)
{
	callTest<TimeSliceSourceTester>( lE, sC, common::sourcetypes::timeslice, iterations);
}

void testPointSource(const LatticeExtents lE, const common::sourcecontents sC, const int iterations)
{
	callTest<PointSourceTester>( lE, sC, common::sourcetypes::point, iterations);
}

void testPsCorrelator(const LatticeExtents lE, const CorrelatorDirection cD, const SpinorFillTypes sF, const ReferenceValues rV)
{
	callTest<ComponentwiseCorrelatorTester>("ps", lE, cD, sF, rV);
}

void testAvpsCorrelator(const LatticeExtents lE, const CorrelatorDirection cD, const SpinorFillTypes sF, const ReferenceValues rV)
{
	callTest<ComponentwiseCorrelatorTester>("avps", lE, cD, sF, rV);
}

void testScCorrelator(const LatticeExtents lE, const CorrelatorDirection cD, const SpinorFillTypes sF, const ReferenceValues rV)
{
	callTest<ColorwiseCorrelatorTester>("sc", lE, cD, sF, rV);
}

void testVxCorrelator(const LatticeExtents lE, const CorrelatorDirection cD, const SpinorFillTypes sF, const ReferenceValues rV)
{
	callTest<ColorwiseCorrelatorTester>("vx", lE, cD, sF, rV);
}

void testVyCorrelator(const LatticeExtents lE, const CorrelatorDirection cD, const SpinorFillTypes sF, const ReferenceValues rV)
{
	callTest<ColorwiseCorrelatorTester>("vy", lE, cD, sF, rV);
}

void testVzCorrelator(const LatticeExtents lE, const CorrelatorDirection cD, const SpinorFillTypes sF, const ReferenceValues rV)
{
	callTest<ColorwiseCorrelatorTester>("vz", lE, cD, sF, rV);
}

void testAxCorrelator(const LatticeExtents lE, const CorrelatorDirection cD, const SpinorFillTypes sF, const ReferenceValues rV)
{
	callTest<ColorwiseCorrelatorTester>("ax", lE, cD, sF, rV);
}

void testAyCorrelator(const LatticeExtents lE, const CorrelatorDirection cD, const SpinorFillTypes sF, const ReferenceValues rV)
{
	callTest<ColorwiseCorrelatorTester>("ay", lE, cD, sF, rV);
}

void testAzCorrelator(const LatticeExtents lE, const CorrelatorDirection cD, const SpinorFillTypes sF, const ReferenceValues rV)
{
	callTest<ColorwiseCorrelatorTester>("az", lE, cD, sF, rV);
}

BOOST_AUTO_TEST_SUITE(SRC_VOLUME)

	BOOST_AUTO_TEST_CASE( SRC_VOLUME_1 )
	{
		testVolumeSource( LatticeExtents{4,4}, common::sourcecontents::one, 12 );
	}

	BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(SRC_VOLUME_2, 2)

	BOOST_AUTO_TEST_CASE( SRC_VOLUME_2 )
	{
		testVolumeSource( LatticeExtents{4,8}, common::sourcecontents::gaussian, 2000);
	}

	BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(SRC_VOLUME_3, 2)

	BOOST_AUTO_TEST_CASE( SRC_VOLUME_3 )
	{
		testVolumeSource( LatticeExtents{8,4}, common::sourcecontents::z4, 1000 );
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SRC_ZSLICE)

	BOOST_AUTO_TEST_CASE( SRC_ZSLICE_1 )
	{
		testZSliceSource( LatticeExtents{4,4}, common::sourcecontents::one, 1 );
	}

	BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(SRC_ZSLICE_2, 2)

	BOOST_AUTO_TEST_CASE( SRC_ZSLICE_2 )
	{
		testZSliceSource( LatticeExtents{4,8}, common::sourcecontents::gaussian, 1000);
	}

	BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(SRC_ZSLICE_3, 2)

	BOOST_AUTO_TEST_CASE( SRC_ZSLICE_3 )
	{
		testZSliceSource( LatticeExtents{16,4}, common::sourcecontents::z4, 20 );
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SRC_TIMESLICE)

	BOOST_AUTO_TEST_CASE( SRC_TIMESLICE_1 )
	{
		testTimeSliceSource( LatticeExtents{4,4}, common::sourcecontents::one, 1 );
	}

	BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(SRC_TIMESLICE_2, 2)

	BOOST_AUTO_TEST_CASE( SRC_TIMESLICE_2 )
	{
		testTimeSliceSource( LatticeExtents{4,8}, common::sourcecontents::gaussian, 2000);
	}

	BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES(SRC_TIMESLICE_3, 2)

	BOOST_AUTO_TEST_CASE( SRC_TIMESLICE_3 )
	{
		testTimeSliceSource( LatticeExtents{16,4}, common::sourcecontents::z4, 20 );
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SRC_POINT)

	BOOST_AUTO_TEST_CASE( SRC_POINT_1 )
	{
		testPointSource( LatticeExtents{4,4}, common::sourcecontents::one, 12 );
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_PS_Z)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testPsCorrelator(LatticeExtents{ ns8, nt4}, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::zero}, ReferenceValues(ns8, 0));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testPsCorrelator(LatticeExtents{ ns8, nt4}, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::one}, ReferenceValues(ns8, 48));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_PS_T)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testPsCorrelator(LatticeExtents {ns8, nt4}, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::zero}, ReferenceValues(nt4, 0));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testPsCorrelator(LatticeExtents {ns8, nt4}, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::one}, ReferenceValues(nt4, 48));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_SC_Z)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testScCorrelator(LatticeExtents{ ns8, nt4}, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(ns8, 0));
	}

	BOOST_AUTO_TEST_CASE( zero2 )
	{
		testScCorrelator(LatticeExtents{ ns8, nt4}, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(ns8, 0));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testScCorrelator(LatticeExtents{ ns8, nt4}, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::ascendingReal, SpinorFillType::oneZero, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(ns8, 1872));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_SC_T)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testScCorrelator(LatticeExtents{ ns4, nt8}, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(nt8, 0));
	}

	BOOST_AUTO_TEST_CASE( zero2 )
	{
		testScCorrelator(LatticeExtents{ ns4, nt8}, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(nt8, 0));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testScCorrelator(LatticeExtents{ ns4, nt8}, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::ascendingReal, SpinorFillType::oneZero, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(nt8, 1872));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_VX_Z)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testVxCorrelator(LatticeExtents{ ns8, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(ns8, 0));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testVxCorrelator(LatticeExtents{ ns8, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(ns8, 192.));
	}

	BOOST_AUTO_TEST_CASE( nonZero2 )
	{
		testVxCorrelator(LatticeExtents{ ns8, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::oneZero, SpinorFillType::zeroOne, SpinorFillType::oneZero, SpinorFillType::zeroOne}, ReferenceValues(ns8, 96.));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_VX_T)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testVxCorrelator(LatticeExtents{ ns8, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(nt12, 0));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testVxCorrelator(LatticeExtents{ ns8, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(nt12, 192.));
	}

	BOOST_AUTO_TEST_CASE( nonZero2 )
	{
		testVxCorrelator(LatticeExtents{ ns8, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::oneZero, SpinorFillType::zeroOne, SpinorFillType::oneZero, SpinorFillType::zeroOne}, ReferenceValues(nt12, 96.));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_VY_Z)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testVyCorrelator(LatticeExtents{ ns8, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(ns8, 0));
	}

	BOOST_AUTO_TEST_CASE( zero2 )
	{
		testVyCorrelator(LatticeExtents{ ns8, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(ns8, 0.));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testVyCorrelator(LatticeExtents{ ns8, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::oneZero, SpinorFillType::zeroOne, SpinorFillType::oneZero, SpinorFillType::zeroOne}, ReferenceValues(ns8, 96.));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_VY_T)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testVyCorrelator(LatticeExtents{ ns8, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(nt12, 0));
	}

	BOOST_AUTO_TEST_CASE( zero2 )
	{
		testVyCorrelator(LatticeExtents{ ns8, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(nt12, 0.));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testVyCorrelator(LatticeExtents{ ns8, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::oneZero, SpinorFillType::zeroOne, SpinorFillType::oneZero, SpinorFillType::zeroOne}, ReferenceValues(nt12, 96.));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_VZ_Z)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testVzCorrelator(LatticeExtents{ ns4, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(ns4, 0));
	}

	BOOST_AUTO_TEST_CASE( zero2 )
	{
		testVzCorrelator(LatticeExtents{ ns4, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(ns4, 0.));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testVzCorrelator(LatticeExtents{ ns4, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::oneZero, SpinorFillType::zeroOne, SpinorFillType::oneZero, SpinorFillType::zeroOne}, ReferenceValues(ns4, 96.));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_VZ_T)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testVzCorrelator(LatticeExtents{ ns12, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(nt12, 0));
	}

	BOOST_AUTO_TEST_CASE( zero2 )
	{
		testVzCorrelator(LatticeExtents{ ns12, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(nt12, 0.));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testVzCorrelator(LatticeExtents{ ns12, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::oneZero, SpinorFillType::zeroOne, SpinorFillType::oneZero, SpinorFillType::zeroOne}, ReferenceValues(nt12, 96.));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_AX_Z)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testAxCorrelator(LatticeExtents{ ns4, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(ns4, 0));
	}

	BOOST_AUTO_TEST_CASE( zero2 )
	{
		testAxCorrelator(LatticeExtents{ ns4, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(ns4, 0.));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testAxCorrelator(LatticeExtents{ ns4, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::ascendingReal, SpinorFillType::oneZero, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(ns4, 144.));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_AX_T)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testAxCorrelator(LatticeExtents{ ns12, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(nt12, 0));
	}

	BOOST_AUTO_TEST_CASE( zero2 )
	{
		testAxCorrelator(LatticeExtents{ ns12, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(nt12, 0.));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testAxCorrelator(LatticeExtents{ ns12, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::ascendingReal, SpinorFillType::oneZero, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(nt12, 144.));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_AY_Z)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testAyCorrelator(LatticeExtents{ ns4, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(ns4, 0));
	}

	BOOST_AUTO_TEST_CASE( zero2 )
	{
		testAyCorrelator(LatticeExtents{ ns4, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(ns4, 0.));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testAyCorrelator(LatticeExtents{ ns4, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::ascendingReal, SpinorFillType::oneZero, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(ns4, -144.));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_AY_T)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testAyCorrelator(LatticeExtents{ ns12, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(nt12, 0));
	}

	BOOST_AUTO_TEST_CASE( zero2 )
	{
		testAyCorrelator(LatticeExtents{ ns12, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(nt12, 0.));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testAyCorrelator(LatticeExtents{ ns12, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::ascendingReal, SpinorFillType::oneZero, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(nt12, -144.));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_AZ_Z)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testAzCorrelator(LatticeExtents{ ns4, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(ns4, 0));
	}

	BOOST_AUTO_TEST_CASE( zero2 )
	{
		testAzCorrelator(LatticeExtents{ ns4, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(ns4, 0.));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testAzCorrelator(LatticeExtents{ ns4, nt12 }, CorrelatorDirection::spatialZ, SpinorFillTypes{SpinorFillType::ascendingReal, SpinorFillType::oneZero, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(ns4, -432.));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_AZ_T)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testAzCorrelator(LatticeExtents{ ns12, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero, SpinorFillType::zero}, ReferenceValues(nt12, 0));
	}

	BOOST_AUTO_TEST_CASE( zero2 )
	{
		testAzCorrelator(LatticeExtents{ ns12, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::one, SpinorFillType::one, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(nt12, 0.));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testAzCorrelator(LatticeExtents{ ns12, nt12 }, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::ascendingReal, SpinorFillType::oneZero, SpinorFillType::one, SpinorFillType::one}, ReferenceValues(nt12, -432.));
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CORRELATOR_AVPS_T)

	BOOST_AUTO_TEST_CASE( zero1 )
	{
		testAvpsCorrelator(LatticeExtents {ns8, nt8}, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::zero}, ReferenceValues(nt8, 0));
	}

	BOOST_AUTO_TEST_CASE( nonZero1 )
	{
		testAvpsCorrelator(LatticeExtents {ns8, nt8}, CorrelatorDirection::temporal, SpinorFillTypes{SpinorFillType::one}, ReferenceValues(nt8, -48));
	}

BOOST_AUTO_TEST_SUITE_END()

