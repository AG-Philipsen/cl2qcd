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

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE OPENCL_MODULE_SPINORS
#include <boost/test/unit_test.hpp>

#include "SpinorTester.hpp"

const ReferenceValues calculateReferenceValues_saxpyAndSquarenormEvenOdd(const int latticeVolume, const hmc_complex alphaIn)
{
	hmc_complex coeff = {1. - alphaIn.re, 0. - alphaIn.im};
	return ReferenceValues { (coeff.im * coeff.im + coeff.re * coeff.re) * latticeVolume * sumOfIntegersSquared(24)};
}

struct SaxpyAndSquarenormEvenOddTestParameters: public SpinorTestParameters
{
	SaxpyAndSquarenormEvenOddTestParameters(const LatticeExtents lE, const hmc_complex cN):
		TestParameters(lE), SpinorTestParameters(lE, SpinorFillTypes{SpinorFillType::ascendingComplex}), coefficient(cN) {};
	const hmc_complex coefficient;
};

struct SaxpyAndSquarenormEvenOddTester: public EvenOddSpinorTester
{
	SaxpyAndSquarenormEvenOddTester(const ParameterCollection & pC, const SaxpyAndSquarenormEvenOddTestParameters & tP):
		EvenOddSpinorTester("saxpy_AND_squarenorm_eo", pC, tP, calculateReferenceValues_saxpyAndSquarenormEvenOdd(calculateEvenOddSpinorfieldSize(tP.latticeExtents), tP.coefficient))
	{
		const hardware::buffers::Spinor in(tP.latticeExtents, device);
		const hardware::buffers::Spinor in2(tP.latticeExtents, device);
		const hardware::buffers::Spinor out(tP.latticeExtents, device);
		const hardware::buffers::Plain<hmc_complex> sqnorm(1, device);
		const hardware::buffers::Plain<hmc_complex> complexNum(1, device);

		EvenOddSpinorfieldCreator sf(tP.latticeExtents);
		in.load(sf.createSpinorfield(tP.fillTypes.at(0)));
		in2.load(sf.createSpinorfield(tP.fillTypes.at(0)));
		complexNum.load(&tP.coefficient);

		code->saxpy_AND_squarenorm_eo_device(&in, &in2, &complexNum, &out, &sqnorm);

		hmc_complex cpu_res = {0., 0.};
		sqnorm.dump(&cpu_res);
		kernelResult[0] = cpu_res.re;
	}
};

void testSaxpyAndSquarenormEvenOdd(const LatticeExtents lE, const hmc_complex cN)
{
	SaxpyAndSquarenormEvenOddTestParameters parametersForThisTest{lE, cN};
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt, true);
	hardware::code::OpenClKernelParametersMockupForMergedSpinorKernels kernelParameters(parametersForThisTest.ns, parametersForThisTest. nt);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	SaxpyAndSquarenormEvenOddTester( parameterCollection, parametersForThisTest);
}

BOOST_AUTO_TEST_SUITE(SF_SAXPY_AND_SQUARENORM_EO)

BOOST_AUTO_TEST_CASE( SF_SAXPY_AND_SQUARENORM_EO_1 )
{
	testSaxpyAndSquarenormEvenOdd( LatticeExtents {ns4,nt4}, hmc_complex{0.,0.});
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_AND_SQUARENORM_EO_2 )
{
	testSaxpyAndSquarenormEvenOdd( LatticeExtents {ns4,nt4}, hmc_complex{1.,0.});
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_AND_SQUARENORM_EO_3 )
{
	testSaxpyAndSquarenormEvenOdd( LatticeExtents {ns4,nt4}, hmc_complex{-1.,0.});
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_AND_SQUARENORM_EO_4 )
{
	testSaxpyAndSquarenormEvenOdd( LatticeExtents {ns4,nt4}, hmc_complex{0.,1.});
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_AND_SQUARENORM_EO_5 )
{
	testSaxpyAndSquarenormEvenOdd( LatticeExtents {ns4,nt4}, hmc_complex{0.,-1.});
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_AND_SQUARENORM_EO_6 )
{
	testSaxpyAndSquarenormEvenOdd( LatticeExtents {ns4,nt4}, hmc_complex{1.,1.});
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_AND_SQUARENORM_EO_7 )
{
	testSaxpyAndSquarenormEvenOdd( LatticeExtents {ns4,nt4}, hmc_complex{-1.,1.});
}

BOOST_AUTO_TEST_CASE( SF_SAXPY_AND_SQUARENORM_EO_8 )
{
	testSaxpyAndSquarenormEvenOdd( LatticeExtents {ns4,nt4}, hmc_complex{-1.,-1.});
}
BOOST_AUTO_TEST_SUITE_END()
