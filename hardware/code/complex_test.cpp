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
		TestParameters(lE), alpha(alphaIn), beta(betaIn), gamma(1.) {}
	hmc_complex alpha, beta;
	hmc_float gamma;
};

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

/**
 * @todo: What is the purpose of "mulitple_operation"?
 * If this is meaningful, it must be covered in the tests!
 */

struct ComplexProductTester: public ComplexTester
{
	ComplexProductTester(const ParameterCollection pC, const ComplexTestParameters tP) :
	   ComplexTester("product", pC, tP, ReferenceValues{nonTrivialParameter, nonTrivialParameter} )
	{
		code->set_complex_to_product_device(alpha, beta, result);
	}
};

struct ComplexRatioTester: public ComplexTester{
	ComplexRatioTester(const ParameterCollection pC, const ComplexTestParameters tP) :
	   ComplexTester("ratio", pC, tP, ReferenceValues{nonTrivialParameter, nonTrivialParameter} )
	{
		code->set_complex_to_ratio_device(alpha, beta, result);
	}
};

struct ComplexSumTester: public ComplexTester{
	ComplexSumTester(const ParameterCollection pC, const ComplexTestParameters tP) :
	   ComplexTester("sum", pC, tP, ReferenceValues{nonTrivialParameter, nonTrivialParameter} )
	{
		code->set_complex_to_sum_device(alpha, beta, result);
	}
};

struct ComplexDifferenceTester: public ComplexTester
{
	ComplexDifferenceTester(const ParameterCollection pC, const ComplexTestParameters tP) :
	   ComplexTester("difference", pC, tP, ReferenceValues{nonTrivialParameter, nonTrivialParameter} )
	{
		code->set_complex_to_difference_device(alpha, beta, result);
	}
};

struct ComplexConvertTester: public ComplexTester
{
	ComplexConvertTester(const ParameterCollection pC, const ComplexTestParameters tP) :
		ComplexTester("convert", pC, tP, ReferenceValues{nonTrivialParameter, nonTrivialParameter} )
	{
		hardware::buffers::Plain<hmc_float> gamma(1, device);
		hmc_float tmp = tP.gamma;
		gamma.load(&tmp);
		code->set_complex_to_float_device(&gamma, result);
	}
};

template<class TesterClass>
void callTest(const LatticeExtents lE)
{
	ComplexTestParameters parametersForThisTest(lE, {1.,0.}, {1.,0.}); //@todo: make adjustable
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt);
	hardware::code::OpenClKernelParametersMockup kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, parametersForThisTest);
}

void testProduct(const LatticeExtents lE)
{
	callTest<ComplexProductTester>(lE);
}

void testRatio(const LatticeExtents lE)
{
	callTest<ComplexRatioTester>(lE);
}

void testSum(const LatticeExtents lE)
{
	callTest<ComplexSumTester>(lE);
}

void testDifference(const LatticeExtents lE)
{
	callTest<ComplexDifferenceTester>(lE);
}

void testConvert(const LatticeExtents lE)
{
	callTest<ComplexConvertTester>(lE);
}

/**
 * @todo: most of the tests should be covered with one test only...
 */

BOOST_AUTO_TEST_SUITE(PRODUCT)

	//@todo: add more tests like in "product_input_{1-12}"
	BOOST_AUTO_TEST_CASE( PRODUCT_1 )
	{
		testProduct(LatticeExtents{ns4, nt4});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(RATIO)

	//@todo: add more tests like in "ratio_input_{1-8}"
	BOOST_AUTO_TEST_CASE( RATIO_1 )
	{
		testRatio(LatticeExtents{ns4, nt4});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SUM)

	//@todo: add more tests like in "sum_input_{1-6}"
	BOOST_AUTO_TEST_CASE( RATIO_1 )
	{
		testSum(LatticeExtents{ns4, nt4});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(DIFFERENCE)

	//@todo: add more tests like in "difference_input_{1-6}"
	BOOST_AUTO_TEST_CASE( RATIO_1 )
	{
		testSum(LatticeExtents{ns4, nt4});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(CONVERT)

	//@todo: add more tests like in "convert_input_1,2"
	BOOST_AUTO_TEST_CASE( CONVERT_1 )
	{
	    testConvert(LatticeExtents{ns4,nt4});
	}
	
BOOST_AUTO_TEST_SUITE_END()
