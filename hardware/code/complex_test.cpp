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
#define BOOST_TEST_MODULE OPENCL_MODULE_COMPLEX

#include "kernelTester.hpp"
#include "complex.hpp"

struct ComplexTestParameters : public TestParameters
{
	ComplexTestParameters(const LatticeExtents lE, const hmc_complex alphaIn, const hmc_complex betaIn):
	TestParameters(lE), alpha(alphaIn), beta(betaIn), gamma(1.) {};
	ComplexTestParameters(const hmc_complex alphaIn, const hmc_complex betaIn):
	TestParameters(LatticeExtents{ns4,nt4}), alpha(alphaIn), beta(betaIn), gamma(1.) {};
	ComplexTestParameters(const hmc_float gammaIn):
	TestParameters(LatticeExtents{ns4,nt4}), alpha({1.,1.}), beta({1.,1.}), gamma(gammaIn) {};
	hmc_complex alpha, beta;
	hmc_float gamma;
};

ReferenceValues calculateReferenceValues_product(const hmc_complex alpha, const hmc_complex beta)
{
	return ReferenceValues{alpha.re * beta.re - alpha.im * beta.im, alpha.im * beta.re + alpha.re * beta.im};
}

ReferenceValues calculateReferenceValues_ratio(const hmc_complex alpha, const hmc_complex beta)
{
	return ReferenceValues{(alpha.re * beta.re + alpha.im * beta.im) / (beta.re * beta.re + beta.im * beta.im), (alpha.im * beta.re - alpha.re * beta.im) / (beta.re * beta.re + beta.im * beta.im)};
}

ReferenceValues calculateReferenceValues_sum(const hmc_complex alpha, const hmc_complex beta)
{
	return ReferenceValues{alpha.re + beta.re, alpha.im + beta.im};
}

ReferenceValues calculateReferenceValues_difference(const hmc_complex alpha, const hmc_complex beta)
{
	return ReferenceValues{alpha.re - beta.re, alpha.im - beta.im};
}

ReferenceValues calculateReferenceValues_convert(const hmc_float gamma)
{
	return ReferenceValues{gamma, 0.};
}

struct ComplexTester : public KernelTester
{
	ComplexTester(std::string kernelName, const ParameterCollection pC, const ComplexTestParameters tP, const ReferenceValues rV) :
	    KernelTester(kernelName, pC.hardwareParameters, pC.kernelParameters, tP, rV)
	{
		code = device->getComplexCode();
		alpha = new hardware::buffers::Plain<hmc_complex>(1, device);
		beta = new hardware::buffers::Plain<hmc_complex>(1, device);
		result = new hardware::buffers::Plain<hmc_complex>(1, device);
		alpha->load(&tP.alpha);
		beta->load(&tP.beta);
	}
	
	void storeResultAsComplex(){
		hmc_complex tmp;
		result->dump(&tmp);
		kernelResult.at(0) = tmp.re;
		kernelResult.at(1) = tmp.im;
	}
	
	virtual ~ComplexTester(){
		storeResultAsComplex();
		delete alpha;
		delete beta;
		delete result;
		code = nullptr;
	}
    
   protected:
	const hardware::code::Complex * code;
	hardware::buffers::Plain<hmc_complex> *alpha;
	hardware::buffers::Plain<hmc_complex> *beta;
	hardware::buffers::Plain<hmc_complex> *result;
};

struct ComplexProductTester: public ComplexTester
{
	ComplexProductTester(const ParameterCollection pC, const ComplexTestParameters tP) :
	ComplexTester("product", pC, tP, calculateReferenceValues_product(tP.alpha, tP.beta))
	{
		code->set_complex_to_product_device(alpha, beta, result);
	}
};

struct ComplexRatioTester: public ComplexTester{
	ComplexRatioTester(const ParameterCollection pC, const ComplexTestParameters tP) :
	ComplexTester("ratio", pC, tP, calculateReferenceValues_ratio(tP.alpha, tP.beta) )
	{
		code->set_complex_to_ratio_device(alpha, beta, result);
	}
};

struct ComplexSumTester: public ComplexTester{
	ComplexSumTester(const ParameterCollection pC, const ComplexTestParameters tP) :
	ComplexTester("sum", pC, tP, calculateReferenceValues_sum(tP.alpha, tP.beta) )
	{
		code->set_complex_to_sum_device(alpha, beta, result);
	}
};

struct ComplexDifferenceTester: public ComplexTester
{
	ComplexDifferenceTester(const ParameterCollection pC, const ComplexTestParameters tP) :
	ComplexTester("difference", pC, tP, calculateReferenceValues_difference(tP.alpha, tP.beta) )
	{
		code->set_complex_to_difference_device(alpha, beta, result);
	}
};

struct ComplexConvertTester: public ComplexTester
{
	ComplexConvertTester(const ParameterCollection pC, const ComplexTestParameters tP) :
	ComplexTester("convert", pC, tP, calculateReferenceValues_convert(tP.gamma) )
	{
		hardware::buffers::Plain<hmc_float> gamma(1, device);
		hmc_float tmp = tP.gamma;
		gamma.load(&tmp);
		code->set_complex_to_float_device(&gamma, result);
	}
};

template<class TesterClass>
void callTest(const hmc_float gamma)
{
	ComplexTestParameters parametersForThisTest(gamma);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt);
	hardware::code::OpenClKernelParametersMockup kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, parametersForThisTest);
}


template<class TesterClass>
void callTest(const hmc_complex alpha, const hmc_complex beta)
{
	ComplexTestParameters parametersForThisTest(alpha, beta);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt);
	hardware::code::OpenClKernelParametersMockup kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, parametersForThisTest);
}

void testProduct(const hmc_complex alpha, const hmc_complex beta)
{
	callTest<ComplexProductTester>(alpha, beta);
}

void testRatio(const hmc_complex alpha, const hmc_complex beta)
{
	callTest<ComplexRatioTester>(alpha, beta);
}

void testSum(const hmc_complex alpha, const hmc_complex beta)
{
	callTest<ComplexSumTester>(alpha, beta);
}

void testDifference(const hmc_complex alpha, const hmc_complex beta)
{
	callTest<ComplexDifferenceTester>(alpha, beta);
}

void testConvert(const hmc_float gamma)
{
	callTest<ComplexConvertTester>(gamma);
}

BOOST_AUTO_TEST_SUITE(PRODUCT)

	BOOST_AUTO_TEST_CASE( PRODUCT_1 )
	{
		testProduct({-nonTrivialParameter, nonTrivialParameter}, {nonTrivialParameter, nonTrivialParameter});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(RATIO)

	BOOST_AUTO_TEST_CASE( RATIO_1 )
	{
		testRatio({nonTrivialParameter, nonTrivialParameter}, {-nonTrivialParameter, nonTrivialParameter});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SUM)

	BOOST_AUTO_TEST_CASE( SUM_1 )
	{
		testSum({nonTrivialParameter, nonTrivialParameter}, {nonTrivialParameter, nonTrivialParameter});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(DIFFERENCE)

	BOOST_AUTO_TEST_CASE( DIFFERENCE_1 )
	{
		testDifference({nonTrivialParameter, nonTrivialParameter}, {nonTrivialParameter, nonTrivialParameter});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CONVERT)

	BOOST_AUTO_TEST_CASE( CONVERT_1 )
	{
		testConvert(nonTrivialParameter);
	}

BOOST_AUTO_TEST_SUITE_END()
