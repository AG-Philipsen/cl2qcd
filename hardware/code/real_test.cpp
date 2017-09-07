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
#define BOOST_TEST_MODULE OPENCL_MODULE_REAL

#include "kernelTester.hpp"
#include "real.hpp"

struct RealTestParameters : public TestParameters
{
	RealTestParameters(const LatticeExtents lE, const hmc_float alpha, const hmc_float beta) :
		TestParameters(lE), alpha(alpha), beta(beta) {};
	RealTestParameters(const hmc_float alpha, const hmc_float beta) :
		TestParameters(LatticeExtents{ns4,nt4}), alpha(alpha), beta(beta) {};
	double alpha, beta;
};

struct TestParametersRealUpdate : public RealTestParameters
{
	double gamma, delta, epsilon, zeta;

	TestParametersRealUpdate(const double a, const double b, const double c, const double d, const double e, const double f = 1):
		RealTestParameters(a, b), gamma(c), delta(d), epsilon(e), zeta(f) {}
};

ReferenceValues calculateReferenceValues_product(const hmc_float alpha, const hmc_float beta)
{
	if (alpha == 0. and beta == 0.)
		return ReferenceValues{0.};
	if (alpha == -nonTrivialParameter and beta == nonTrivialParameter)
		return ReferenceValues{-1.*nonTrivialParameter * nonTrivialParameter};
	return defaultReferenceValues();
}

ReferenceValues calculateReferenceValues_ratio(const hmc_float alpha, const hmc_float beta)
{
	if (alpha == 2.*nonTrivialParameter and beta == -nonTrivialParameter)
		return ReferenceValues{-2.};
	return defaultReferenceValues();
}

ReferenceValues calculateReferenceValues_sum(const hmc_float alpha, const hmc_float beta)
{
	if (alpha == nonTrivialParameter and beta == -nonTrivialParameter)
		return ReferenceValues{0.};
	if (alpha == 3.*nonTrivialParameter and beta == -5.*nonTrivialParameter)
		return ReferenceValues{-2.*nonTrivialParameter};
	return defaultReferenceValues();
}

ReferenceValues calculateReferenceValues_difference(const hmc_float alpha, const hmc_float beta)
{
	if (alpha == nonTrivialParameter and beta == nonTrivialParameter)
		return ReferenceValues{0.};
	if (alpha == 3.*nonTrivialParameter and beta == -5.*nonTrivialParameter)
		return ReferenceValues{8.*nonTrivialParameter};
	return defaultReferenceValues();
}

ReferenceValues calculateReferenceValues_realUpdate(TestParametersRealUpdate tP)
{
	return ReferenceValues{-tP.alpha * tP.beta * tP.gamma / (tP.delta * tP.epsilon)};
}

ReferenceValues calculateReferenceValues_realUpdateBeta(TestParametersRealUpdate tP)
{
	return ReferenceValues{-tP.alpha * tP.beta / tP.gamma};
}

ReferenceValues calculateReferenceValues_realUpdateZeta(TestParametersRealUpdate tP)
{
	return ReferenceValues{(tP.alpha * tP.beta * tP.gamma) / (tP.delta * tP.epsilon * (tP.beta - tP.alpha) + tP.beta * tP.gamma * (1. - tP.zeta * tP.delta))};
}

struct RealTester : public KernelTester
{
	RealTester(std::string kernelName, const ParameterCollection pC, const RealTestParameters tP, const ReferenceValues rV) :
		KernelTester(kernelName, pC.hardwareParameters, pC.kernelParameters, tP, rV)
	{
		code = device->getRealCode();

		alpha = new hardware::buffers::Plain<hmc_float>(1, device);
		beta = new hardware::buffers::Plain<hmc_float>(1, device);
		result = new hardware::buffers::Plain<hmc_float>(1, device);
		alpha->load(&tP.alpha);
		beta->load(&tP.beta);
	}

	virtual ~RealTester(){
		result->dump(&kernelResult.at(0));
	}

   protected:
	const hardware::code::Real * code;
	hardware::buffers::Plain<hmc_float> *alpha;
	hardware::buffers::Plain<hmc_float> *beta;
	hardware::buffers::Plain<hmc_float> *result;
};

struct RealProductTester: public RealTester
{
	RealProductTester(const ParameterCollection pC, const RealTestParameters tP) :
			RealTester("product",pC, tP, calculateReferenceValues_product(tP.alpha, tP.beta))
	{
		code->set_real_to_product_device(alpha, beta, result);
	}
};

struct RealRatioTester: public RealTester
{
	RealRatioTester(const ParameterCollection pC, const RealTestParameters tP) :
	   RealTester("ratio", pC, tP, calculateReferenceValues_ratio(tP.alpha, tP.beta) )
	{
		code->set_real_to_ratio_device(alpha, beta, result);
	}
};

struct RealSumTester: public RealTester
{
	RealSumTester(const ParameterCollection pC, const RealTestParameters tP) :
	   RealTester("sum", pC, tP, calculateReferenceValues_sum(tP.alpha, tP.beta) ){
		code->set_real_to_sum_device(alpha, beta, result);
	}
};

struct RealDifferenceTester: public RealTester
{
	RealDifferenceTester(const ParameterCollection pC, const RealTestParameters tP) :
	   RealTester("difference",pC, tP, calculateReferenceValues_difference(tP.alpha, tP.beta) )
	{
		code->set_real_to_difference_device(alpha, beta, result);
	}
};

struct RealSetVectorElementTester: public KernelTester
{
	RealSetVectorElementTester(const ParameterCollection pC, const TestParameters tP) :
		KernelTester("Access_vector_element", pC.hardwareParameters, pC.kernelParameters, tP, ReferenceValues{0.})
	{
		std::vector<hmc_float> vector_host(4);

		const hardware::code::Real * code = device->getRealCode();
		hardware::buffers::Plain<hmc_float> scalar_buf(1, device);
		hardware::buffers::Plain<hmc_float> vector_buf(4, device);

		for(uint i = 0; i<vector_host.size(); i++)
		{
			vector_host.at(i) = i*nonTrivialParameter;
		}
		vector_buf.load(&vector_host.at(0));
		for(uint i=0; i<vector_host.size(); i++)
		{
			code->set_real_to_vector_element_device(&vector_buf, i, &scalar_buf);
			hmc_float cpu_res;
			scalar_buf.dump(&cpu_res);
			BOOST_REQUIRE_CLOSE(cpu_res, vector_host.at(i), 1.e-8);
		}
	}
};

struct RealAccessVectorElementTester: public KernelTester
{
	RealAccessVectorElementTester(const ParameterCollection pC, const TestParameters tP) :
		KernelTester("Access_vector_element", pC.hardwareParameters, pC.kernelParameters, tP, ReferenceValues{0.})
	{
		std::vector<hmc_float> vector_host(4);

		const hardware::code::Real * code = device->getRealCode();
		hardware::buffers::Plain<hmc_float> scalar_buf(1, device);
		hardware::buffers::Plain<hmc_float> vector_buf(4, device);
		hmc_float scalar_host;

	   scalar_host = nonTrivialParameter;
	   scalar_buf.load(&scalar_host);
	   for(uint i=0; i<vector_host.size(); i++){
		  code->set_vector_element_to_real_device(&scalar_buf, i, &vector_buf);
		  std::vector<hmc_float> cpu_res(4);
		  vector_buf.dump(&cpu_res.at(0));
		  BOOST_REQUIRE_CLOSE(cpu_res[i], scalar_host, 1.e-8);
	   }
	}
};

struct RealUpdateTesterAlpha: public RealTester
{
	RealUpdateTesterAlpha(const ParameterCollection pC, const TestParametersRealUpdate tP) :
		RealTester("update", pC, tP, calculateReferenceValues_realUpdate(tP) )
	{
	   hardware::buffers::Plain<hmc_float> gamma(1, device);
	   hardware::buffers::Plain<hmc_float> delta(1, device);
	   hardware::buffers::Plain<hmc_float> epsilon(1, device);

	   gamma.load(&tP.gamma);
	   delta.load(&tP.delta);
	   epsilon.load(&tP.epsilon);

	   code->update_alpha_cgm_device(alpha, beta, &gamma, &delta, &epsilon, 1, result);
	}
};

struct RealUpdateTesterBeta: public RealTester
{
	RealUpdateTesterBeta(const ParameterCollection pC, const TestParametersRealUpdate tP) :
		RealTester("update", pC, tP, calculateReferenceValues_realUpdateBeta(tP) )
	{
	   hardware::buffers::Plain<hmc_float> gamma(1, device);
	   gamma.load(&tP.gamma);
	   code->update_beta_cgm_device(alpha, beta, &gamma, 1, result);
	}
};

struct RealUpdateTesterZeta: public RealTester
{
	RealUpdateTesterZeta(const ParameterCollection pC, const TestParametersRealUpdate tP) :
		RealTester("update", pC, tP, calculateReferenceValues_realUpdateZeta(tP) )
	{
	   hardware::buffers::Plain<hmc_float> gamma(1, device);
	   hardware::buffers::Plain<hmc_float> delta(1, device);
	   hardware::buffers::Plain<hmc_float> epsilon(1, device);
	   hardware::buffers::Plain<hmc_float> zeta(1, device);

	   gamma.load(&tP.gamma);
	   delta.load(&tP.delta);
	   epsilon.load(&tP.epsilon);
	   zeta.load(&tP.zeta);

	   code->update_zeta_cgm_device(alpha, beta, &gamma, &delta, &epsilon, &zeta, 1, result);
	}
};

template<class TesterClass>
void callTest(const hmc_float alpha, const hmc_float beta)
{
	RealTestParameters parametersForThisTest(alpha, beta);
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.latticeExtents);
	hardware::code::OpenClKernelParametersMockup kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, parametersForThisTest);
}

void testProduct(const hmc_float alpha, const hmc_float beta)
{
	callTest<RealProductTester>(alpha, beta);
}

void testRatio(const hmc_float alpha, const hmc_float beta)
{
	callTest<RealRatioTester>(alpha, beta);
}

void testSum(const hmc_float alpha, const hmc_float beta)
{
	callTest<RealSumTester>(alpha, beta);
}

void testDifference(const hmc_float alpha, const hmc_float beta)
{
	callTest<RealDifferenceTester>(alpha, beta);
}

template<class TesterClass>
void testAccess()
{
	TestParameters parametersForThisTest(LatticeExtents{ns4, nt4});
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt);
	hardware::code::OpenClKernelParametersMockup kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, parametersForThisTest);
}

template<class TesterClass>
void testUpdate(const std::vector<hmc_float> valuesIn)
{
	TestParametersRealUpdate parametersForThisTest(valuesIn.at(0), valuesIn.at(1), valuesIn.at(2), valuesIn.at(3), valuesIn.at(4), valuesIn.at(5));
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns, parametersForThisTest.nt);
	hardware::code::OpenClKernelParametersMockup kernelParameters(parametersForThisTest.ns, parametersForThisTest.nt);
	ParameterCollection parameterCollection{hardwareParameters, kernelParameters};
	TesterClass(parameterCollection, parametersForThisTest);
}

BOOST_AUTO_TEST_SUITE(PRODUCT)

	BOOST_AUTO_TEST_CASE( PRODUCT_1 )
	{
		testProduct(0., 0.);
	}

	BOOST_AUTO_TEST_CASE( PRODUCT_2 )
	{
		testProduct(-nonTrivialParameter, nonTrivialParameter);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(RATIO)

	BOOST_AUTO_TEST_CASE( RATIO_1 )
	{
		testRatio(2.*nonTrivialParameter, -nonTrivialParameter);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SUM)

	BOOST_AUTO_TEST_CASE( SUM_1 )
	{
		testSum(nonTrivialParameter, -nonTrivialParameter);
	}

	BOOST_AUTO_TEST_CASE( SUM_2 )
	{
		testSum(3.*nonTrivialParameter, -5.*nonTrivialParameter);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(DIFFERENCE)

	BOOST_AUTO_TEST_CASE( DIFFERENCE_1 )
	{
		testDifference(nonTrivialParameter, nonTrivialParameter);
	}

	BOOST_AUTO_TEST_CASE( DIFFERENCE_2 )
	{
		testDifference(3.*nonTrivialParameter, -5.*nonTrivialParameter);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(SET_ELEMENT)

	BOOST_AUTO_TEST_CASE( SET_ELEMENT_1 )
	{
		testAccess<RealSetVectorElementTester>();
	}

	BOOST_AUTO_TEST_CASE( SET_ELEMENT_2 )
	{
		testAccess<RealSetVectorElementTester>();
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(ACCESS_ELEMENT)

	BOOST_AUTO_TEST_CASE( ACCESS_ELEMENT_1 )
	{
		testAccess<RealAccessVectorElementTester>();
	}

	BOOST_AUTO_TEST_CASE( ACCESS_ELEMENT_2 )
	{
		testAccess<RealAccessVectorElementTester>();
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(REAL_UPDATE)

	BOOST_AUTO_TEST_CASE( ALPHA_1 )
	{
		testUpdate<RealUpdateTesterAlpha>(std::vector<hmc_float>{ 1.35, -1.25, 3.4, -2., 4.5, 0.});
	}

	BOOST_AUTO_TEST_CASE( BETA_1 )
	{
		testUpdate<RealUpdateTesterBeta>(std::vector<hmc_float>{ 1.35, -1.25, 3.4, 0., 0., 0.});
	}

	BOOST_AUTO_TEST_CASE( ZETA_1 )
	{
		testUpdate<RealUpdateTesterZeta>(std::vector<hmc_float>{ 1.35, -1.25, 3.4, -2., 4.5, 2.3});
	}

BOOST_AUTO_TEST_SUITE_END()

