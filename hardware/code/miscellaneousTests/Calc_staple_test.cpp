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
#define BOOST_TEST_MODULE staple_test
#include <boost/test/unit_test.hpp>

#include "../GaugefieldTester.hpp"

struct TestCode : public hardware::code::Opencl_Module
{
	virtual size_t get_read_write_size(const std::string&) const
	{
		return 0;
	};
	virtual uint64_t get_flop_size(const std::string&) const
	{
		return 0;
	};

	TestCode(const hardware::code::OpenClKernelParametersInterface & kP, hardware::Device * device) :
		Opencl_Module(kP, device)
	{
		testKernel = createKernel("staple_test") << get_device()->getGaugefieldCode()->get_sources() << "../hardware/code/miscellaneousTests/staple_test.cl";
	};
	virtual ~TestCode()
	{
		clReleaseKernel(testKernel);
	};

	void runTestKernel(const hardware::buffers::SU3 * gf, const hardware::buffers::Plain<hmc_float> * out, const int gs, const int ls)
	{
		err = clSetKernelArg(testKernel, 0, sizeof(cl_mem), gf->get_cl_buffer());
		BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
		err = clSetKernelArg(testKernel, 1, sizeof(cl_mem), out->get_cl_buffer());
		BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);

		get_device()->enqueue_kernel(testKernel, gs, ls);
	}
	cl_kernel testKernel;
	cl_int err;
};

struct StapleTester : public GaugefieldTester
{
	StapleTester(const ParameterCollection pC, const GaugefieldTestParameters tP, const ReferenceValues rV):
		GaugefieldTester("StapleTest", pC, tP, rV)
	{
		testCode = new TestCode(pC.kernelParameters, device);
		hardware::buffers::Plain<hmc_float> out (calculateGaugefieldSize(tP.latticeExtents), device);

		if(device->get_device_type() == CL_DEVICE_TYPE_GPU)
		{
			gs = calculateGaugefieldSize(tP.latticeExtents);
			ls = 64;
		} else {
			gs = device->get_num_compute_units();
			ls = 1;
		}

		testCode->runTestKernel(GaugefieldTester::gaugefieldBuffer, &out, gs, ls);

		host_out = new hmc_float[calculateGaugefieldSize(tP.latticeExtents)];
		out.dump(host_out);

		hmc_float result = 0;
		for(int i = 0; i < calculateGaugefieldSize(tP.latticeExtents); i++) {
			result += host_out[i];
		}
		kernelResult.at(0) = result;
	}
	~StapleTester()
	{
		if(testCode)
			delete testCode;
		if(host_out)
			delete host_out;
	}
	TestCode * testCode;
	int ls, gs;
	hmc_float * host_out;
};

BOOST_AUTO_TEST_CASE( STAPLE_TEST )
{
	GaugefieldTestParameters parametersForThisTest {LatticeExtents{4,4}, GaugefieldFillType::cold};
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns,parametersForThisTest.nt);
	hardware::code::OpenClKernelParametersMockup kernelParameters(parametersForThisTest.ns,parametersForThisTest.nt);
	ParameterCollection parameterCollection(hardwareParameters, kernelParameters);
	StapleTester tester(parameterCollection, parametersForThisTest, defaultReferenceValues());
}
