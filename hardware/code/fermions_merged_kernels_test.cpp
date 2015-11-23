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
	else if ( coefficient.re == -nonTrivialMassParameter and coefficient.im == 2.*nonTrivialMassParameter)
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
	if (massParameters.kappa == -4.*nonTrivialMassParameter and massParameters.mu == 11*nonTrivialMassParameter)
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
	if (massParameters.kappa == 2.*nonTrivialMassParameter and massParameters.mu == -5*nonTrivialMassParameter)
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
		if (massParameters.kappa == 3*nonTrivialMassParameter and massParameters.mu == -7*nonTrivialMassParameter)
		{
			return ReferenceValues{ latticeVolume * (-630.5134172436968) };
		}
	}
	else if (gF == GaugefieldFillType::nonTrivial and massParameters.kappa == 3*nonTrivialMassParameter and massParameters.mu == -7*nonTrivialMassParameter)
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
		if (massParameters.kappa == -2*nonTrivialMassParameter and massParameters.mu == -8*nonTrivialMassParameter)
		{
			return ReferenceValues{ latticeVolume * (478.7145794216685) };
		}
	}
	else if (gF == GaugefieldFillType::nonTrivial and massParameters.kappa == -2*nonTrivialMassParameter and massParameters.mu == -8*nonTrivialMassParameter)
	{
		return ReferenceValues{ latticeVolume * (191.461710119754) };
	}
	return defaultReferenceValues();
}

struct SaxpyAndGamma5EvenOddTestParameters: public EvenOddFermionTestParameters
{
	SaxpyAndGamma5EvenOddTestParameters(const LatticeExtents latticeExtentsIn, const SpinorFillTypes spinorFillTypesIn, const hmc_complex coefficientIn) :
		EvenOddFermionTestParameters(calculateReferenceValues_saxpyAndGamma5EvenOdd(calculateEvenOddSpinorfieldSize(latticeExtentsIn), coefficientIn),
				latticeExtentsIn, spinorFillTypesIn), coefficient(coefficientIn) {};
	const hmc_complex coefficient;
};

struct SaxpyAndGamma5EvenOddTester : public FermionTester
{
	SaxpyAndGamma5EvenOddTester(const ParameterCollection & pc, const SaxpyAndGamma5EvenOddTestParameters & testParameters) :
		FermionTester("saxpy_AND_gamma5", pc, testParameters)
	{
		const hardware::buffers::Spinor in(spinorfieldEvenOddElements, device);
		const hardware::buffers::Spinor in2(spinorfieldEvenOddElements, device);
		const hardware::buffers::Spinor out(spinorfieldEvenOddElements, device);

		in.load(createSpinorfield(testParameters.fillTypes.at(0)));
		in2.load(createSpinorfield(testParameters.fillTypes.at(1)));
		out.load(createSpinorfield(SpinorFillType::zero));

		code->saxpy_AND_gamma5_eo_device(&in, &in2, testParameters.coefficient, &out);

		spinor * sf_in;
		sf_in = new spinor[spinorfieldEvenOddElements];
		out.dump(sf_in);
		kernelResult[0] = count_sf(sf_in, spinorfieldEvenOddElements);
		delete sf_in;
	}
};

/**
 * @todo: move the counting part to the destructor of the base class
 */
struct MTmSitediagonalMinusAndGamma5EvenOddTester : public FermionTester
{
	MTmSitediagonalMinusAndGamma5EvenOddTester(const ParameterCollection & pc, const EvenOddTwistedMassTestParameters & testParameters) :
		FermionTester("m_tm_sitediagonal_minus_AND_gamma5_eo", pc, testParameters)
	{
		const hardware::buffers::Spinor in(spinorfieldEvenOddElements, device);
		const hardware::buffers::Spinor out(spinorfieldEvenOddElements, device);

		in.load(createSpinorfield(testParameters.fillTypes.at(0)));
		out.load(createSpinorfield(SpinorFillType::zero));

		code->M_tm_sitediagonal_minus_AND_gamma5_eo_device(&in, &out, testParameters.massParameters.getMubar());

		spinor * sf_in;
		sf_in = new spinor[spinorfieldEvenOddElements];
		out.dump(sf_in);
		kernelResult[0] = count_sf(sf_in, spinorfieldEvenOddElements);
		delete sf_in;
	}
};

struct MTmSitediagonalAndGamma5EvenOddTester : public FermionTester
{
	MTmSitediagonalAndGamma5EvenOddTester(const ParameterCollection & pc, const EvenOddTwistedMassTestParameters & testParameters) :
		FermionTester("m_tm_sitediagonal_AND_gamma5_eo", pc, testParameters)
	{
		const hardware::buffers::Spinor in(spinorfieldEvenOddElements, device);
		const hardware::buffers::Spinor out(spinorfieldEvenOddElements, device);

		in.load(createSpinorfield(testParameters.fillTypes.at(0)));
		out.load(createSpinorfield(SpinorFillType::zero));

		code->M_tm_sitediagonal_AND_gamma5_eo_device(&in, &out, testParameters.massParameters.getMubar());

		spinor * sf_in;
		sf_in = new spinor[spinorfieldEvenOddElements];
		out.dump(sf_in);
		kernelResult[0] = count_sf(sf_in, spinorfieldEvenOddElements);
		delete sf_in;
	}
};

struct DslashEvenOddAndMTmInverseSitediagonalMinusTester : public FermionmatrixEvenOddTester
{
	DslashEvenOddAndMTmInverseSitediagonalMinusTester(const ParameterCollection & pc, const EvenOddTwistedMassTestParameters & testParameters, const bool evenOrOdd) :
		FermionmatrixEvenOddTester("dslash_AND_m_tm_inverse_sitediagonal_minus", pc, testParameters)
	{
		in->load(createSpinorfield(testParameters.fillTypes.at(0)));
		out->load(createSpinorfield(SpinorFillType::zero));

		code->dslash_AND_M_tm_inverse_sitediagonal_minus_eo_device(in, out, gaugefieldBuffer, evenOrOdd, testParameters.massParameters.kappa,
				testParameters.massParameters.getMubar());

		spinor * sf_in;
		sf_in = new spinor[spinorfieldEvenOddElements];
		out->dump(sf_in);
		kernelResult[0] = count_sf(sf_in, spinorfieldEvenOddElements);
		delete sf_in;
	}
};

struct DslashEvenOddAndMTmInverseSitediagonalTester : public FermionmatrixEvenOddTester
{
	DslashEvenOddAndMTmInverseSitediagonalTester(const ParameterCollection & pc, const EvenOddTwistedMassTestParameters & testParameters, const bool evenOrOdd) :
		FermionmatrixEvenOddTester("dslash_AND_m_tm_inverse_sitediagonal", pc, testParameters)
	{
		in->load(createSpinorfield(testParameters.fillTypes.at(0)));
		out->load(createSpinorfield(SpinorFillType::zero));

		code->dslash_AND_M_tm_inverse_sitediagonal_eo_device(in, out, gaugefieldBuffer, evenOrOdd, testParameters.massParameters.kappa,
				testParameters.massParameters.getMubar());

		spinor * sf_in;
		sf_in = new spinor[spinorfieldEvenOddElements];
		out->dump(sf_in);
		kernelResult[0] = count_sf(sf_in, spinorfieldEvenOddElements);
		delete sf_in;
	}
};

/**
 * @todo: combine these calls using templates
 */
void testSaxpyAndGamma5EvenOdd(const LatticeExtents lE, const SpinorFillTypes sF, const hmc_complex coefficient)
{
	SaxpyAndGamma5EvenOddTestParameters parametersForThisTest{lE, sF, coefficient};
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.SpinorTestParameters::ns, parametersForThisTest.SpinorTestParameters::nt, parametersForThisTest.needEvenOdd);
	hardware::code::OpenClKernelParametersMockupForMergedFermionKernels kernelParameters(parametersForThisTest.SpinorTestParameters::ns, parametersForThisTest.SpinorTestParameters::nt);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	SaxpyAndGamma5EvenOddTester( parameterCollection, parametersForThisTest);
}

void testMTmSitediagonalMinusAndGamma5EvenOdd(const LatticeExtents lE, const SpinorFillType sF, const TwistedMassMassParameters mP)
{
	EvenOddTwistedMassTestParameters parametersForThisTest{lE, sF, mP, calculateReferenceValues_mTmSitediagonalMinusAndGamma5EvenOdd};
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.SpinorTestParameters::ns, parametersForThisTest.SpinorTestParameters::nt, parametersForThisTest.needEvenOdd);
	hardware::code::OpenClKernelParametersMockupForTwistedMass kernelParameters(parametersForThisTest.SpinorTestParameters::ns,
			parametersForThisTest.SpinorTestParameters::nt, true, true);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	MTmSitediagonalMinusAndGamma5EvenOddTester( parameterCollection, parametersForThisTest);
}

void testMTmSitediagonalAndGamma5EvenOdd(const LatticeExtents lE, const SpinorFillType sF, const TwistedMassMassParameters mP)
{
	EvenOddTwistedMassTestParameters parametersForThisTest{lE, sF, mP, calculateReferenceValues_mTmSitediagonalAndGamma5EvenOdd};
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.SpinorTestParameters::ns, parametersForThisTest.SpinorTestParameters::nt, parametersForThisTest.needEvenOdd);
	hardware::code::OpenClKernelParametersMockupForTwistedMass kernelParameters(parametersForThisTest.SpinorTestParameters::ns,
			parametersForThisTest.SpinorTestParameters::nt, true, true);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	MTmSitediagonalAndGamma5EvenOddTester( parameterCollection, parametersForThisTest);
}

void testDslashEvenOddAndMTmInverseSitediagonalMinus(const LatticeExtents lE, const SpinorFillType sF, const GaugefieldFillType gF, const TwistedMassMassParameters mP, const bool evenOrOdd)
{
	EvenOddTwistedMassTestParameters parametersForThisTest{lE, sF, gF, mP, calculateReferenceValues_dslashEvenOddAndMTmInverseSitediagonalMinus};
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.SpinorTestParameters::ns, parametersForThisTest.SpinorTestParameters::nt, parametersForThisTest.needEvenOdd);
	hardware::code::OpenClKernelParametersMockupForTwistedMass kernelParameters(parametersForThisTest.SpinorTestParameters::ns,
			parametersForThisTest.SpinorTestParameters::nt, true, true);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	DslashEvenOddAndMTmInverseSitediagonalMinusTester( parameterCollection, parametersForThisTest, evenOrOdd);
}

void testDslashEvenOddAndMTmInverseSitediagonal(const LatticeExtents lE, const SpinorFillType sF, const GaugefieldFillType gF, const TwistedMassMassParameters mP, const bool evenOrOdd)
{
	EvenOddTwistedMassTestParameters parametersForThisTest{lE, sF, gF, mP, calculateReferenceValues_dslashEvenOddAndMTmInverseSitediagonal};
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.SpinorTestParameters::ns, parametersForThisTest.SpinorTestParameters::nt, parametersForThisTest.needEvenOdd);
	hardware::code::OpenClKernelParametersMockupForTwistedMass kernelParameters(parametersForThisTest.SpinorTestParameters::ns,
			parametersForThisTest.SpinorTestParameters::nt, true, true);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	DslashEvenOddAndMTmInverseSitediagonalTester( parameterCollection, parametersForThisTest, evenOrOdd);
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
		testDslashEvenOddAndMTmInverseSitediagonal(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, TwistedMassMassParameters{-2*nonTrivialMassParameter, -8*nonTrivialMassParameter}, true);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_EO_4)
	{
		testDslashEvenOddAndMTmInverseSitediagonal(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, TwistedMassMassParameters{-2*nonTrivialMassParameter, -8*nonTrivialMassParameter}, false);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_EO_5)
	{
		testDslashEvenOddAndMTmInverseSitediagonal(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::nonTrivial, TwistedMassMassParameters{-2*nonTrivialMassParameter, -8*nonTrivialMassParameter}, true);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_EO_6)
	{
		testDslashEvenOddAndMTmInverseSitediagonal(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::nonTrivial, TwistedMassMassParameters{-2*nonTrivialMassParameter, -8*nonTrivialMassParameter}, false);
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
		testDslashEvenOddAndMTmInverseSitediagonalMinus(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, TwistedMassMassParameters{3.*nonTrivialMassParameter, -7*nonTrivialMassParameter}, true);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_MINUS_EO_4)
	{
		testDslashEvenOddAndMTmInverseSitediagonalMinus(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, TwistedMassMassParameters{3.*nonTrivialMassParameter, -7*nonTrivialMassParameter}, false);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_MINUS_EO_5)
	{
		testDslashEvenOddAndMTmInverseSitediagonalMinus(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::nonTrivial, TwistedMassMassParameters{3.*nonTrivialMassParameter, -7*nonTrivialMassParameter}, true);
	}

	BOOST_AUTO_TEST_CASE( DSLASH_AND_M_TM_INVERSE_SITEDIAGONAL_MINUS_EO_6)
	{
		testDslashEvenOddAndMTmInverseSitediagonalMinus(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::nonTrivial, TwistedMassMassParameters{3.*nonTrivialMassParameter, -7*nonTrivialMassParameter}, false);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(M_TM_SITEDIAGONAL_AND_GAMMA5_EO )

	BOOST_AUTO_TEST_CASE(M_TM_SITEDIAGONAL_AND_GAMMA5_EO_1)
	{
		testMTmSitediagonalAndGamma5EvenOdd(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, TwistedMassMassParameters{0.,0.});
	}

	BOOST_AUTO_TEST_CASE(M_TM_SITEDIAGONAL_AND_GAMMA5_EO_2)
	{
		testMTmSitediagonalAndGamma5EvenOdd(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, TwistedMassMassParameters{-4.*nonTrivialMassParameter,11.*nonTrivialMassParameter});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(M_TM_SITEDIAGONAL_MINUS_AND_GAMMA5_EO )

	BOOST_AUTO_TEST_CASE(M_TM_SITEDIAGONAL_MINUS_AND_GAMMA5_EO_1)
	{
		testMTmSitediagonalMinusAndGamma5EvenOdd(LatticeExtents{ns4, nt8}, SpinorFillType::ascendingComplex, TwistedMassMassParameters{0.,0.});
	}

	BOOST_AUTO_TEST_CASE(M_TM_SITEDIAGONAL_MINUS_AND_GAMMA5_EO_2)
	{
		testMTmSitediagonalMinusAndGamma5EvenOdd(LatticeExtents{ns4, nt8}, SpinorFillType::ascendingComplex, TwistedMassMassParameters{2.*nonTrivialMassParameter,-5.*nonTrivialMassParameter});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SAXPY_AND_GAMMA5_EO )

	BOOST_AUTO_TEST_CASE(SAXPY_AND_GAMMA5_EO_1)
	{
		testSaxpyAndGamma5EvenOdd(LatticeExtents{ns4, nt8}, SpinorFillTypes{SpinorFillType::ascendingComplex, SpinorFillType::ascendingComplex}, hmc_complex{0.,0.});
	}

	BOOST_AUTO_TEST_CASE(SAXPY_AND_GAMMA5_EO_2)
	{
		testSaxpyAndGamma5EvenOdd(LatticeExtents{ns4, nt12}, SpinorFillTypes{SpinorFillType::ascendingComplex, SpinorFillType::ascendingComplex}, hmc_complex{-nonTrivialMassParameter, 2.*nonTrivialMassParameter});
	}
	
BOOST_AUTO_TEST_SUITE_END()
