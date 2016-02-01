/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
 *   2015 Christopher Pinke
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
#define BOOST_TEST_MODULE OPENCL_MODULE_FERMIONS_MERGED_KERNELS
#include <boost/test/unit_test.hpp>

#include "FermionTester.hpp"

const ReferenceValues calculateReferenceValues_saxpyAndGamma5EvenOdd(const int latticeVolume, const hmc_complex coefficient)
{
	if( coefficient.re == 0. and coefficient.im == 0.)
	{
		return ReferenceValues{ latticeVolume * (-144.) };
	}
	else if ( coefficient.re == -nonTrivialParameter and coefficient.im == 2.*nonTrivialParameter)
	{
		return ReferenceValues{ latticeVolume * (-161.7776640000001) };
	}
	return defaultReferenceValues();
}

const ReferenceValues calculateReferenceValues_mTmSitediagonalAndGamma5EvenOdd(const int latticeVolume, const TwistedMassMassParameters massParameters)
{
	if (massParameters.kappa == 0. and massParameters.mu == 0.)
	{
		return ReferenceValues{ latticeVolume * (-144.) };
	}
	if (massParameters.kappa == -4.*nonTrivialParameter and massParameters.mu == 11*nonTrivialParameter)
	{
		return ReferenceValues{ latticeVolume * (-127.905098563584) };
	}
	return defaultReferenceValues();
}

const ReferenceValues calculateReferenceValues_mTmSitediagonalMinusAndGamma5EvenOdd(const int latticeVolume, const TwistedMassMassParameters massParameters)
{
	if (massParameters.kappa == 0. and massParameters.mu == 0.)
	{
		return ReferenceValues{ latticeVolume * (-144.) };
	}
	if (massParameters.kappa == 2.*nonTrivialParameter and massParameters.mu == -5*nonTrivialParameter)
	{
		return ReferenceValues{ latticeVolume * (-147.65793214464) };
	}
	return defaultReferenceValues();
}

const ReferenceValues calculateReferenceValues_dslashEvenOddAndMTmInverseSitediagonalMinus(const int latticeVolume, const GaugefieldFillType gF, const TwistedMassMassParameters massParameters)
{
	if (gF == GaugefieldFillType::cold )
	{
		if (massParameters.kappa == 1. and massParameters.mu == 0.)
		{
			return ReferenceValues{ latticeVolume * (-2400.) };
		}
		if (massParameters.kappa == 3*nonTrivialParameter and massParameters.mu == -7*nonTrivialParameter)
		{
			return ReferenceValues{ latticeVolume * (-630.5134172436968) };
		}
	}
	else if (gF == GaugefieldFillType::nonTrivial and massParameters.kappa == 3*nonTrivialParameter and massParameters.mu == -7*nonTrivialParameter)
	{
		return ReferenceValues{ latticeVolume * (-262.5581454299608) };
	}
	return defaultReferenceValues();
}

const ReferenceValues calculateReferenceValues_dslashEvenOddAndMTmInverseSitediagonal(const int latticeVolume, const GaugefieldFillType gF, const TwistedMassMassParameters massParameters)
{
	if (gF == GaugefieldFillType::cold )
	{
		if (massParameters.kappa == 1. and massParameters.mu == 0.)
		{
			return ReferenceValues{ latticeVolume * (-2400.) };
		}
		if (massParameters.kappa == -2*nonTrivialParameter and massParameters.mu == -8*nonTrivialParameter)
		{
			return ReferenceValues{ latticeVolume * (478.7145794216685) };
		}
	}
	else if (gF == GaugefieldFillType::nonTrivial and massParameters.kappa == -2*nonTrivialParameter and massParameters.mu == -8*nonTrivialParameter)
	{
		return ReferenceValues{ latticeVolume * (191.461710119754) };
	}
	return defaultReferenceValues();
}

struct SaxpyAndGamma5EvenOddTestParameters: public FermionTestParameters
{
	SaxpyAndGamma5EvenOddTestParameters(const LatticeExtents latticeExtentsIn, const SpinorFillTypes spinorFillTypesIn, const hmc_complex coefficientIn) :
		TestParameters(latticeExtentsIn), FermionTestParameters(latticeExtentsIn, spinorFillTypesIn), coefficient(coefficientIn) {};
	const hmc_complex coefficient;
};

struct SaxpyAndGamma5EvenOddTester : public FermionmatrixTesterWithSumAsKernelResult<EvenOddFermionmatrixTester>
{
	SaxpyAndGamma5EvenOddTester(const ParameterCollection & pc, const SaxpyAndGamma5EvenOddTestParameters & tP) :
		FermionmatrixTesterWithSumAsKernelResult<EvenOddFermionmatrixTester>("saxpy_AND_gamma5", pc, tP, calculateReferenceValues_saxpyAndGamma5EvenOdd(calculateEvenOddSpinorfieldSize(tP.SpinorTestParameters::latticeExtents), tP.coefficient))
	{
		EvenOddSpinorfieldCreator sf(tP.SpinorTestParameters::latticeExtents);
		const hardware::buffers::Spinor in2(tP.SpinorTestParameters::latticeExtents, device);
		in2.load(sf.createSpinorfield(tP.fillTypes.at(1)));
		code->saxpy_AND_gamma5_eo_device(in, &in2, tP.coefficient, out);
	}
};

struct MTmSitediagonalMinusAndGamma5EvenOddTester : public FermionmatrixTesterWithSumAsKernelResult<EvenOddFermionmatrixTester>
{
	MTmSitediagonalMinusAndGamma5EvenOddTester(const ParameterCollection & pc, const TwistedMassTestParameters & tP) :
		FermionmatrixTesterWithSumAsKernelResult<EvenOddFermionmatrixTester>("m_tm_sitediagonal_minus_AND_gamma5_eo", pc, tP,
				calculateReferenceValues_mTmSitediagonalMinusAndGamma5EvenOdd(calculateEvenOddSpinorfieldSize(tP.SpinorTestParameters::latticeExtents), tP.massParameters) )
	{
		code->M_tm_sitediagonal_minus_AND_gamma5_eo_device(in, out, tP.massParameters.getMubar());
	}
};

struct MTmSitediagonalAndGamma5EvenOddTester : public FermionmatrixTesterWithSumAsKernelResult<EvenOddFermionmatrixTester>
{
	MTmSitediagonalAndGamma5EvenOddTester(const ParameterCollection & pc, const TwistedMassTestParameters & tP) :
		FermionmatrixTesterWithSumAsKernelResult<EvenOddFermionmatrixTester>("m_tm_sitediagonal_AND_gamma5_eo", pc, tP,
				calculateReferenceValues_mTmSitediagonalAndGamma5EvenOdd(calculateEvenOddSpinorfieldSize(tP.SpinorTestParameters::latticeExtents), tP.massParameters))
	{
		code->M_tm_sitediagonal_AND_gamma5_eo_device(in, out, tP.massParameters.getMubar());
	}
};

struct DslashEvenOddAndMTmInverseSitediagonalMinusTester : public FermionmatrixTesterWithSumAsKernelResult<EvenOddFermionmatrixTester>
{
	DslashEvenOddAndMTmInverseSitediagonalMinusTester(const ParameterCollection & pc, const TwistedMassTestParameters & tP, const bool evenOrOdd) :
		FermionmatrixTesterWithSumAsKernelResult<EvenOddFermionmatrixTester>("dslash_AND_m_tm_inverse_sitediagonal_minus", pc, tP,
				calculateReferenceValues_dslashEvenOddAndMTmInverseSitediagonalMinus(calculateEvenOddSpinorfieldSize(tP.SpinorTestParameters::latticeExtents), tP.GaugefieldTestParameters::fillType, tP.massParameters) )
	{
		code->dslash_AND_M_tm_inverse_sitediagonal_minus_eo_device(in, out, gaugefieldBuffer, evenOrOdd, tP.massParameters.kappa,
				tP.massParameters.getMubar());
	}
};

struct DslashEvenOddAndMTmInverseSitediagonalTester : public FermionmatrixTesterWithSumAsKernelResult<EvenOddFermionmatrixTester>
{
	DslashEvenOddAndMTmInverseSitediagonalTester(const ParameterCollection & pc, const TwistedMassTestParameters & tP, const bool evenOrOdd) :
		FermionmatrixTesterWithSumAsKernelResult<EvenOddFermionmatrixTester>("dslash_AND_m_tm_inverse_sitediagonal", pc, tP,
				calculateReferenceValues_dslashEvenOddAndMTmInverseSitediagonal(calculateEvenOddSpinorfieldSize(tP.SpinorTestParameters::latticeExtents), tP.GaugefieldTestParameters::fillType, tP.massParameters) )
	{
		code->dslash_AND_M_tm_inverse_sitediagonal_eo_device(in, out, gaugefieldBuffer, evenOrOdd, tP.massParameters.kappa,
				tP.massParameters.getMubar());
	}
};

void testSaxpyAndGamma5EvenOdd(const LatticeExtents lE, const SpinorFillTypes sF, const hmc_complex coefficient)
{
	SaxpyAndGamma5EvenOddTestParameters parametersForThisTest{lE, sF, coefficient};
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt, true);
	hardware::code::OpenClKernelParametersMockupForMergedFermionKernels kernelParameters(parametersForThisTest.ns, parametersForThisTest. nt);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	SaxpyAndGamma5EvenOddTester( parameterCollection, parametersForThisTest);
}

template <class TesterClass>
void runTest(const LatticeExtents lE, const SpinorFillType sF, const TwistedMassMassParameters mP)
{
	TwistedMassTestParameters parametersForThisTest{lE, sF, mP};
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt, true);
	hardware::code::OpenClKernelParametersMockupForTwistedMass kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt, true, true);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass( parameterCollection, parametersForThisTest);
}

template <class TesterClass>
void runTest(const LatticeExtents lE, const SpinorFillType sF, const GaugefieldFillType gF, const TwistedMassMassParameters mP, const bool evenOrOdd)
{
	TwistedMassTestParameters parametersForThisTest{lE, sF, gF, mP};
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt, true);
	hardware::code::OpenClKernelParametersMockupForTwistedMass kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt, true, true);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass( parameterCollection, parametersForThisTest, evenOrOdd);
}

void testMTmSitediagonalMinusAndGamma5EvenOdd(const LatticeExtents lE, const SpinorFillType sF, const TwistedMassMassParameters mP)
{
	runTest<MTmSitediagonalMinusAndGamma5EvenOddTester>( lE, sF, mP );
}

void testMTmSitediagonalAndGamma5EvenOdd(const LatticeExtents lE, const SpinorFillType sF, const TwistedMassMassParameters mP)
{
	runTest<MTmSitediagonalAndGamma5EvenOddTester>( lE, sF, mP );
}

void testDslashEvenOddAndMTmInverseSitediagonalMinus(const LatticeExtents lE, const SpinorFillType sF, const GaugefieldFillType gF, const TwistedMassMassParameters mP, const bool evenOrOdd)
{
	runTest<DslashEvenOddAndMTmInverseSitediagonalMinusTester>( lE, sF, gF, mP, evenOrOdd );
}

void testDslashEvenOddAndMTmInverseSitediagonal(const LatticeExtents lE, const SpinorFillType sF, const GaugefieldFillType gF, const TwistedMassMassParameters mP, const bool evenOrOdd)
{
	runTest<DslashEvenOddAndMTmInverseSitediagonalTester>( lE, sF, gF, mP, evenOrOdd );
}

/**
 * @todo: Add more parameters to dslash-tests (thetaT, thetaS, chem. pot., etc.)
 */
BOOST_AUTO_TEST_SUITE(DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_EO )

	BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_EO_1)
	{
		testDslashEvenOddAndMTmInverseSitediagonal(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, TwistedMassMassParameters{1, 0.}, true);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_EO_2)
	{
		testDslashEvenOddAndMTmInverseSitediagonal(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, TwistedMassMassParameters{1, 0.}, false);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_EO_3)
	{
		testDslashEvenOddAndMTmInverseSitediagonal(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, TwistedMassMassParameters{-2*nonTrivialParameter, -8*nonTrivialParameter}, true);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_EO_4)
	{
		testDslashEvenOddAndMTmInverseSitediagonal(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, TwistedMassMassParameters{-2*nonTrivialParameter, -8*nonTrivialParameter}, false);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_EO_5)
	{
		testDslashEvenOddAndMTmInverseSitediagonal(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::nonTrivial, TwistedMassMassParameters{-2*nonTrivialParameter, -8*nonTrivialParameter}, true);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_EO_6)
	{
		testDslashEvenOddAndMTmInverseSitediagonal(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::nonTrivial, TwistedMassMassParameters{-2*nonTrivialParameter, -8*nonTrivialParameter}, false);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_MINUS_EO )

	BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_MINUS_EO_1)
	{
		testDslashEvenOddAndMTmInverseSitediagonalMinus(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, TwistedMassMassParameters{1, 0.}, true);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_MINUS_EO_2)
	{
		testDslashEvenOddAndMTmInverseSitediagonalMinus(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, TwistedMassMassParameters{1, 0.}, false);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_MINUS_EO_3)
	{
		testDslashEvenOddAndMTmInverseSitediagonalMinus(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, TwistedMassMassParameters{3.*nonTrivialParameter, -7*nonTrivialParameter}, true);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_MINUS_EO_4)
	{
		testDslashEvenOddAndMTmInverseSitediagonalMinus(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, TwistedMassMassParameters{3.*nonTrivialParameter, -7*nonTrivialParameter}, false);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_MINUS_EO_5)
	{
		testDslashEvenOddAndMTmInverseSitediagonalMinus(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::nonTrivial, TwistedMassMassParameters{3.*nonTrivialParameter, -7*nonTrivialParameter}, true);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_MINUS_EO_6)
	{
		testDslashEvenOddAndMTmInverseSitediagonalMinus(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::nonTrivial, TwistedMassMassParameters{3.*nonTrivialParameter, -7*nonTrivialParameter}, false);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(M_TM_SITEDIAGONAL_AND_GAMMA5_EO )

	BOOST_AUTO_TEST_CASE(M_TM_SITEDIAGONAL_AND_GAMMA5_EO_1)
	{
		testMTmSitediagonalAndGamma5EvenOdd(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, TwistedMassMassParameters{0.,0.});
	}

	BOOST_AUTO_TEST_CASE(M_TM_SITEDIAGONAL_AND_GAMMA5_EO_2)
	{
		testMTmSitediagonalAndGamma5EvenOdd(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, TwistedMassMassParameters{-4.*nonTrivialParameter,11.*nonTrivialParameter});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(M_TM_SITEDIAGONAL_MINUS_AND_GAMMA5_EO )

	BOOST_AUTO_TEST_CASE(M_TM_SITEDIAGONAL_MINUS_AND_GAMMA5_EO_1)
	{
		testMTmSitediagonalMinusAndGamma5EvenOdd(LatticeExtents{ns4, nt8}, SpinorFillType::ascendingComplex, TwistedMassMassParameters{0.,0.});
	}

	BOOST_AUTO_TEST_CASE(M_TM_SITEDIAGONAL_MINUS_AND_GAMMA5_EO_2)
	{
		testMTmSitediagonalMinusAndGamma5EvenOdd(LatticeExtents{ns4, nt8}, SpinorFillType::ascendingComplex, TwistedMassMassParameters{2.*nonTrivialParameter,-5.*nonTrivialParameter});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAXPY_AND_GAMMA5_EO )

	BOOST_AUTO_TEST_CASE(SAXPY_AND_GAMMA5_EO_1)
	{
		testSaxpyAndGamma5EvenOdd(LatticeExtents{ns4, nt8}, SpinorFillTypes{SpinorFillType::ascendingComplex, SpinorFillType::ascendingComplex}, hmc_complex{0.,0.});
	}

	BOOST_AUTO_TEST_CASE(SAXPY_AND_GAMMA5_EO_2)
	{
		testSaxpyAndGamma5EvenOdd(LatticeExtents{ns4, nt12}, SpinorFillTypes{SpinorFillType::ascendingComplex, SpinorFillType::ascendingComplex}, hmc_complex{-nonTrivialParameter, 2.*nonTrivialParameter});
	}
	
BOOST_AUTO_TEST_SUITE_END()
