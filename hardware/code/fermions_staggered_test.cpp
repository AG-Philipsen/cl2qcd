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
#define BOOST_TEST_MODULE OPENCL_MODULE_FERMIONS_STAGGERED

#include "SpinorStaggeredTester.hpp"
#include "fermions_staggered.hpp"
#include "FermionTester.hpp"

const ReferenceValues calculateReferenceValues_mStaggered(const LatticeExtents lE, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn, const WilsonMassParameters massParametersIn, const ThetaParameters thetaParametersIn)
{
	int latticeVolume = calculateSpinorfieldSize(lE);
	if (massParametersIn.kappa == nonTrivialParameter )
	{
		if (spinorFillTypeIn == SpinorFillType::ascendingComplex and gaugefieldFillTypeIn == GaugefieldFillType::cold )
		{
			return ReferenceValues{latticeVolume * (183.386965938176 - 136.5 * cos(2 * M_PI * thetaParametersIn.thetaS/lE.getNs()) - 45.5 * cos(2 * M_PI * thetaParametersIn.thetaT/lE.getNt()))};
		}
		else if (spinorFillTypeIn == SpinorFillType::one and gaugefieldFillTypeIn == GaugefieldFillType::nonTrivial )
		{
			return ReferenceValues{latticeVolume * (6.04573031577098 + 1.1986620875976 * cos(2 * M_PI * thetaParametersIn.thetaS/lE.getNs()) + 0.3995540291992 * cos( 2 * M_PI * thetaParametersIn.thetaT/lE.getNt()) + 3.62339720366385 * sin(2 * M_PI * thetaParametersIn.thetaS/lE.getNs()) + 1.20779906788795 * sin( 2 * M_PI * thetaParametersIn.thetaT/lE.getNt()))};
		}
		else if (spinorFillTypeIn == SpinorFillType::ascendingComplex and gaugefieldFillTypeIn == GaugefieldFillType::nonTrivial)
		{
			return ReferenceValues{latticeVolume * (183.3871307172304 + 55.58807312524925 * cos(2 * M_PI * thetaParametersIn.thetaS/lE.getNs()) + 18.52935770841643 * cos( 2 * M_PI * thetaParametersIn.thetaT/lE.getNt()) + 107.8954221508676 * sin(2 * M_PI * thetaParametersIn.thetaS/lE.getNs()) + 35.96514071695587 * sin(2 * M_PI * thetaParametersIn.thetaT/lE.getNt()))};
		}
	}
	return defaultReferenceValues();
}

const ReferenceValues calculateReferenceValues_DksEo(const LatticeExtents lE, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn, const ThetaParameters thetaParametersIn, const ImaginaryChemicalPotential chemPotIn)
{
	int latticeVolume = calculateEvenOddSpinorfieldSize(lE);
	if (spinorFillTypeIn == SpinorFillType::ascendingComplex and gaugefieldFillTypeIn == GaugefieldFillType::cold )
	{
		return ReferenceValues{latticeVolume * (182. - 136.5 * cos(2 * M_PI * thetaParametersIn.thetaS/lE.getNs()) - 45.5 * cos(2 * (chemPotIn + M_PI * thetaParametersIn.thetaT/lE.getNt())))};
	}
	else if (spinorFillTypeIn == SpinorFillType::one and gaugefieldFillTypeIn == GaugefieldFillType::nonTrivial )
	{
		return ReferenceValues{latticeVolume * (6.000006163962978 + 1.1986620875976 * cos(2 * M_PI * thetaParametersIn.thetaS/lE.getNs()) + 0.3995540291992 * cos(2 * (chemPotIn + M_PI * thetaParametersIn.thetaT/lE.getNt())) + 3.62339720366385 * sin(2 * M_PI * thetaParametersIn.thetaS/lE.getNs()) + 1.20779906788795 * sin( 2 * (chemPotIn + M_PI * thetaParametersIn.thetaT/lE.getNt())))};
	}
	else if (spinorFillTypeIn == SpinorFillType::ascendingComplex and gaugefieldFillTypeIn == GaugefieldFillType::nonTrivial)
	{
		return ReferenceValues{latticeVolume * (182.0001647790544 + 55.58807312524925 * cos(2 * M_PI * thetaParametersIn.thetaS/lE.getNs()) + 18.52935770841643 * cos(2 * (chemPotIn + M_PI * thetaParametersIn.thetaT/lE.getNt())) + 107.8954221508676 * sin(2 * M_PI * thetaParametersIn.thetaS/lE.getNs()) + 35.96514071695587 * sin(2 * (chemPotIn + M_PI * thetaParametersIn.thetaT/lE.getNt())))};
	}
	return defaultReferenceValues();
}

struct NonEvenOddStaggeredFermionmatrixTester : public StaggeredFermionmatrixTester<hardware::buffers::Plain<su3vec>, NonEvenOddSpinorStaggeredfieldCreator>
{
	NonEvenOddStaggeredFermionmatrixTester(const std::string kernelName, const ParameterCollection pC, const StaggeredFermionsTestParameters tP, const ReferenceValues rV) :
	    StaggeredFermionmatrixTester<hardware::buffers::Plain<su3vec>, NonEvenOddSpinorStaggeredfieldCreator>(kernelName, pC, tP, rV){}
};

struct MTester : public StaggeredFermionmatrixTesterWithSumAsKernelResult<NonEvenOddStaggeredFermionmatrixTester>
{
	MTester(const ParameterCollection pC, const StaggeredFermionsTestParameters tP) :
		StaggeredFermionmatrixTesterWithSumAsKernelResult<NonEvenOddStaggeredFermionmatrixTester>("M_staggered", pC, tP, calculateReferenceValues_mStaggered(tP.latticeExtents, tP.spinorFillType, tP.gaugefieldFillType, tP.massParameters, tP.thetaParameters))
	{
			code->M_staggered_device(in, out,  gaugefieldBuffer, tP.massParameters.kappa);

	}
};

struct EvenOddStaggeredFermionmatrixTester : public StaggeredFermionmatrixTester<hardware::buffers::SU3vec, EvenOddSpinorStaggeredfieldCreator>
{
	EvenOddStaggeredFermionmatrixTester(std::string kernelName, const ParameterCollection pC, const StaggeredFermionsTestParameters tP, const ReferenceValues rV) :
		StaggeredFermionmatrixTester<hardware::buffers::SU3vec, EvenOddSpinorStaggeredfieldCreator>(kernelName, pC, tP, rV){}
};

struct DksEvenOddTester : public StaggeredFermionmatrixTesterWithSumAsKernelResult<EvenOddStaggeredFermionmatrixTester>
{
	DksEvenOddTester(const ParameterCollection pC, const EvenOddStaggeredFermionsTestParameters tP) :
		StaggeredFermionmatrixTesterWithSumAsKernelResult<EvenOddStaggeredFermionmatrixTester>("DKS_eo", pC, tP, calculateReferenceValues_DksEo(tP.latticeExtents, tP.spinorFillType, tP.gaugefieldFillType, tP.thetaParameters, tP.chemicalPotential))
	{
			tP.evenOrOdd ? code->D_KS_eo_device(in, out,  gaugefieldBuffer, EVEN)
			  : code->D_KS_eo_device(in, out,  gaugefieldBuffer, ODD);
	}
};

template<class TesterClass>
void callTest(const LatticeExtents lE, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn, const WilsonMassParameters massParametersIn, const ThetaParameters thetaIn = {0.,0.} , const bool useEvenOddIn = false)
{
	StaggeredFermionsTestParameters parametersForThisTest(lE, spinorFillTypeIn, gaugefieldFillTypeIn, massParametersIn, thetaIn);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.latticeExtents, useEvenOddIn);
	hardware::code::OpenClKernelParametersMockupForFermionsStaggeredTests kernelParameters(parametersForThisTest.latticeExtents, massParametersIn.kappa, useEvenOddIn, thetaIn.thetaT, thetaIn.thetaS);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, parametersForThisTest);
}

template<class TesterClass>
void callTest(const LatticeExtents lE, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn, const bool useEvenOddIn, const bool evenOrOdd, const ThetaParameters thetaIn = {0.,0.}, const ImaginaryChemicalPotential chemPotIn = 0.)
{
	EvenOddStaggeredFermionsTestParameters parametersForThisTest(lE, spinorFillTypeIn, gaugefieldFillTypeIn, evenOrOdd, thetaIn, chemPotIn);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.latticeExtents, useEvenOddIn);
	hardware::code::OpenClKernelParametersMockupForFermionsStaggeredTests kernelParameters(parametersForThisTest.latticeExtents, useEvenOddIn, thetaIn.thetaT, thetaIn.thetaS, chemPotIn);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, parametersForThisTest);
}

void testMStaggered(const LatticeExtents lE, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn, const WilsonMassParameters massParametersIn)
{
	callTest<MTester>(lE, spinorFillTypeIn, gaugefieldFillTypeIn, massParametersIn);
}

void testMStaggeredWithSpecificBC(const LatticeExtents lE, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn, const WilsonMassParameters massParametersIn, const ThetaParameters thetaIn)
{
	callTest<MTester>(lE, spinorFillTypeIn, gaugefieldFillTypeIn, massParametersIn, thetaIn);

}

void testDksEo(const LatticeExtents lE, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn, const bool evenOrOdd)
{
	callTest<DksEvenOddTester>(lE, spinorFillTypeIn, gaugefieldFillTypeIn, true, evenOrOdd);
}

void testDksEoWithSpecificBC(const LatticeExtents lE, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn, const bool useEvenOddIn, const ThetaParameters thetaIn, const bool evenOrOdd)
{
	callTest<DksEvenOddTester>(lE, spinorFillTypeIn, gaugefieldFillTypeIn, useEvenOddIn, evenOrOdd, thetaIn);
}

void testDksEoWithSpecificBCAndChemicalPotential(const LatticeExtents lE, const SpinorFillType spinorFillTypeIn, const GaugefieldFillType gaugefieldFillTypeIn, const bool useEvenOddIn, const ThetaParameters thetaIn, const ImaginaryChemicalPotential chemPotIn, const bool evenOrOdd)
{
	callTest<DksEvenOddTester>(lE, spinorFillTypeIn, gaugefieldFillTypeIn, useEvenOddIn, evenOrOdd, thetaIn, chemPotIn);
}

BOOST_AUTO_TEST_SUITE( M_MATRIX )

	BOOST_AUTO_TEST_CASE( M_MATRIX_1)
	{
		testMStaggered(LatticeExtents{ns12, nt4}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, WilsonMassParameters{nonTrivialParameter});
	}
	BOOST_AUTO_TEST_CASE( M_MATRIX_2)
	{
		testMStaggered(LatticeExtents{ns4, nt8}, SpinorFillType::one, GaugefieldFillType::nonTrivial, WilsonMassParameters{nonTrivialParameter});
	}
	BOOST_AUTO_TEST_CASE( M_MATRIX_3)
	{
		testMStaggered(LatticeExtents{ns8, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::nonTrivial, WilsonMassParameters{nonTrivialParameter});
	}

	BOOST_AUTO_TEST_CASE( M_MATRIX_BC_1)
	{
		testMStaggeredWithSpecificBC(LatticeExtents{ns4, nt8}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, WilsonMassParameters{nonTrivialParameter}, ThetaParameters{0., nonTrivialParameter});
	}
	BOOST_AUTO_TEST_CASE( M_MATRIX_BC_2)
	{
		testMStaggeredWithSpecificBC(LatticeExtents{ns4, nt4}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, WilsonMassParameters{nonTrivialParameter}, ThetaParameters{nonTrivialParameter, 0.});
	}
	BOOST_AUTO_TEST_CASE( M_MATRIX_BC_3)
	{
		testMStaggeredWithSpecificBC(LatticeExtents{ns8, nt8}, SpinorFillType::one, GaugefieldFillType::nonTrivial, WilsonMassParameters{nonTrivialParameter}, ThetaParameters{nonTrivialParameter, nonTrivialParameter});
	}
	BOOST_AUTO_TEST_CASE( M_MATRIX_BC_4)
	{
		testMStaggeredWithSpecificBC(LatticeExtents{ns4, nt4}, SpinorFillType::one, GaugefieldFillType::nonTrivial, WilsonMassParameters{nonTrivialParameter}, ThetaParameters{0., nonTrivialParameter});
	}
	BOOST_AUTO_TEST_CASE( M_MATRIX_BC_5)
	{
		testMStaggeredWithSpecificBC(LatticeExtents{ns4, nt4}, SpinorFillType::ascendingComplex, GaugefieldFillType::nonTrivial, WilsonMassParameters{nonTrivialParameter}, ThetaParameters{nonTrivialParameter, 0.});
	}
	BOOST_AUTO_TEST_CASE( M_MATRIX_BC_6)
	{
		testMStaggeredWithSpecificBC(LatticeExtents{ns4, nt4}, SpinorFillType::ascendingComplex, GaugefieldFillType::nonTrivial, WilsonMassParameters{nonTrivialParameter}, ThetaParameters{nonTrivialParameter, nonTrivialParameter});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE( DKS_EO )

	BOOST_AUTO_TEST_CASE( DKS_EO_1)
	{
		testDksEo(LatticeExtents{ns4, nt4}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, EVEN);
	}
	BOOST_AUTO_TEST_CASE( DKS_EO_2)
	{
		testDksEo(LatticeExtents{ns4, nt4}, SpinorFillType::one, GaugefieldFillType::nonTrivial, ODD);
	}
	BOOST_AUTO_TEST_CASE( DKS_EO_3)
	{
		testDksEo(LatticeExtents{ns4, nt4}, SpinorFillType::ascendingComplex, GaugefieldFillType::nonTrivial, ODD);
	}
	BOOST_AUTO_TEST_CASE( DKS_EO_BC_1)
	{
		testDksEoWithSpecificBC(LatticeExtents{ns4, nt4}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, true, ThetaParameters{1., nonTrivialParameter}, EVEN);
	}
	BOOST_AUTO_TEST_CASE( DKS_EO_BC_2)
	{
		testDksEoWithSpecificBC(LatticeExtents{ns4, nt4}, SpinorFillType::one, GaugefieldFillType::nonTrivial, true, ThetaParameters{nonTrivialParameter, 1.}, ODD);
	}
	BOOST_AUTO_TEST_CASE( DKS_EO_BC_3)
	{
		testDksEoWithSpecificBC(LatticeExtents{ns4, nt4}, SpinorFillType::ascendingComplex, GaugefieldFillType::nonTrivial, true, ThetaParameters{nonTrivialParameter, nonTrivialParameter}, EVEN);
	}
	BOOST_AUTO_TEST_CASE( DKS_EO_BC_CP_1)
	{
		testDksEoWithSpecificBCAndChemicalPotential(LatticeExtents{ns4, nt4}, SpinorFillType::ascendingComplex, GaugefieldFillType::cold, true, ThetaParameters{1., nonTrivialParameter}, ImaginaryChemicalPotential{nonTrivialParameter}, ODD);
	}
	BOOST_AUTO_TEST_CASE( DKS_EO_BC_CP_2)
	{
		testDksEoWithSpecificBCAndChemicalPotential(LatticeExtents{ns4, nt4}, SpinorFillType::one, GaugefieldFillType::nonTrivial, true, ThetaParameters{nonTrivialParameter, 1.}, ImaginaryChemicalPotential{nonTrivialParameter}, EVEN);
	}
	BOOST_AUTO_TEST_CASE( DKS_EO_BC_CP_3)
	{
		testDksEoWithSpecificBCAndChemicalPotential(LatticeExtents{ns4, nt4}, SpinorFillType::ascendingComplex, GaugefieldFillType::nonTrivial, true, ThetaParameters{nonTrivialParameter, nonTrivialParameter}, ImaginaryChemicalPotential{nonTrivialParameter}, ODD);
	}
	
BOOST_AUTO_TEST_SUITE_END()

