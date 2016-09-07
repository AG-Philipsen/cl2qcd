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
#define BOOST_TEST_MODULE localQ_test
#include <boost/test/unit_test.hpp>

#include "testCode.hpp"

const ReferenceValues calculateReferenceValue_localQ(const LatticeExtents lE)
{
	return ReferenceValues{ 72.00012210960028 * lE.getLatticeVolume() };
}

struct LocalQTestCode : public TestCode
{
	LocalQTestCode(const hardware::code::OpenClKernelParametersInterface & kP, hardware::Device * device):
		TestCode(kP, device)
	{
		testKernel = createKernel("localQ_test") << get_device()->getGaugefieldCode()->get_sources()  << "../hardware/code/miscellaneousTests/localQ_test.cl";
	}

	void runTestKernel(const hardware::buffers::SU3 * gf, const hardware::buffers::Plain<hmc_float> * out, const int gs, const int ls) override
	{
		err = clSetKernelArg(testKernel, 0, sizeof(cl_mem), gf->get_cl_buffer());
		BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
		err = clSetKernelArg(testKernel, 1, sizeof(cl_mem), out->get_cl_buffer());
		BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);

		get_device()->enqueue_kernel(testKernel, gs, ls);
	}
};

struct LocalQTester : public OtherKernelTester
{
	LocalQTester(const ParameterCollection pC, const GaugefieldTestParameters tP):
		OtherKernelTester("local_Q_test", pC, tP, calculateReferenceValue_localQ(tP.latticeExtents))
	{
		testCode = new LocalQTestCode(pC.kernelParameters, device);
		testCode->runTestKernel(OtherKernelTester::gaugefieldBuffer, out, gs, ls);
	}
};

BOOST_AUTO_TEST_CASE( LOCAL_Q )
{
	GaugefieldTestParameters parametersForThisTest {LatticeExtents{4,4}, GaugefieldFillType::nonTrivial};
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns,parametersForThisTest.nt);
	hardware::code::OpenClKernelParametersMockup kernelParameters(parametersForThisTest.ns,parametersForThisTest.nt);
	ParameterCollection parameterCollection(hardwareParameters, kernelParameters);
	LocalQTester tester(parameterCollection, parametersForThisTest);
}
