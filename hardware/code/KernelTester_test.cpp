/*
 * Copyright 2012, 2013, 2014 Lars Zeidlewicz, Christopher Pinke,
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

#include <iostream>
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE KernelTester_test
#include <boost/test/unit_test.hpp>

#include "kernelTester.hpp"

struct DoubleKernelTester : public KernelTester
{
	DoubleKernelTester(std::string kernelNameIn, const hardware::HardwareParametersInterface& hPI,
			const hardware::code::OpenClKernelParametersInterface& kP, struct TestParameters tP, const ReferenceValues rV):
		KernelTester(kernelNameIn, hPI, kP, tP, rV)
	{
		kernelResult.at(0) = 1;
	}
};

struct ComplexKernelTester : public KernelTester
{
	ComplexKernelTester(std::string kernelNameIn, const hardware::HardwareParametersInterface& hPI,
			const hardware::code::OpenClKernelParametersInterface& kP, struct TestParameters tP, const ReferenceValues rV) :
		KernelTester(kernelNameIn, hPI, kP, tP, rV)
	{
		kernelResult.at(0) = 1.;
		kernelResult.at(1) = 2.;
	}
};

BOOST_AUTO_TEST_SUITE ( BUILD )

	BOOST_AUTO_TEST_CASE( BUILD_1 )
	{
		const hardware::HardwareParametersMockup params(4,4);
		const hardware::code::OpenClKernelParametersMockup kernelParameters(4,4,true);
		const struct TestParameters testParams(LatticeExtents(4,4));
		BOOST_CHECK_NO_THROW(KernelTester kernelTester("testMockup", params, kernelParameters, testParams, ReferenceValues{0} ) );
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( DOUBLE )

	BOOST_AUTO_TEST_CASE( TRIVIALKERNEL )
	{
		const hardware::HardwareParametersMockup params(4,4);
		const hardware::code::OpenClKernelParametersMockup kernelParameters(4,4,true);
		const struct TestParameters testParams(LatticeExtents(4,4));

		DoubleKernelTester kernelTester("test", params, kernelParameters, testParams, ReferenceValues{1} );
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( DOUBLE_SMALL )

BOOST_AUTO_TEST_CASE_EXPECTED_FAILURES (FAILING_COMPARISON, 1)

	BOOST_AUTO_TEST_CASE( FAILING_COMPARISON )
	{
		const hardware::HardwareParametersMockup params(4,4);
		const hardware::code::OpenClKernelParametersMockup kernelParameters(4,4,true);
		const struct TestParameters testParams(LatticeExtents(4,4));

		DoubleKernelTester kernelTester("test", params, kernelParameters, testParams, ReferenceValues{0.});
	}

	BOOST_AUTO_TEST_CASE( SUCCEEDING_COMPARISON )
	{
		const hardware::HardwareParametersMockup params(4,4);
		const hardware::code::OpenClKernelParametersMockup kernelParameters(4,4,true);
		const struct TestParameters testParams(LatticeExtents(4,4));

		DoubleKernelTester kernelTester("test", params, kernelParameters, testParams, ReferenceValues{1});
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE ( COMPLEX )

	BOOST_AUTO_TEST_CASE( COMPARE_COMPLEX_NUMBERS )
	{
		const hardware::HardwareParametersMockup params(4,4);
		const hardware::code::OpenClKernelParametersMockup kernelParameters(4,4,true);
		const struct TestParameters testParams(LatticeExtents(4,4));

		ComplexKernelTester kernelTester("test", params, kernelParameters, testParams, ReferenceValues{1.,2.});
	}

BOOST_AUTO_TEST_SUITE_END()
