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
#define BOOST_TEST_MODULE SU2 SU3 Extend
#include <boost/test/unit_test.hpp>

#include "testCode.hpp"

#define NUM_ELEMENTS 1024
#define LOCAL_SIZE 128

struct ExtendTestCode : public TestCode
{
	ExtendTestCode(const hardware::code::OpenClKernelParametersInterface & kP, hardware::Device * device):
		TestCode(kP, device)
	{
		testKernel = createKernel("extendKernel") << get_device()->getGaugefieldCode()->get_sources() << "../hardware/code/miscellaneousTests/su2su3extend.cl";
	}

	void runExtendKernel(const hardware::buffers::Plain<Matrixsu3> * out, const hardware::buffers::Plain<Matrixsu2> * in, const hardware::buffers::Plain<cl_int> * d_rand, cl_ulong elems)
	{
		cl_int err;
		err = clSetKernelArg(testKernel, 0, sizeof(cl_mem), out->get_cl_buffer());
		BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
		err = clSetKernelArg(testKernel, 1, sizeof(cl_mem), in->get_cl_buffer());
		BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
		err = clSetKernelArg(testKernel, 2, sizeof(cl_mem), d_rand->get_cl_buffer());
		BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
		err = clSetKernelArg(testKernel, 3, sizeof(cl_ulong), &elems);
		BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
		get_device()->enqueueKernel(testKernel, NUM_ELEMENTS);
	}
};

struct ExtendTester : public GaugefieldTester
{
	ExtendTester(const ParameterCollection pC, const GaugefieldTestParameters tP):
		GaugefieldTester("su2su3Extend", pC, tP, ReferenceValues{0.})
	{
		for(auto device: system->get_devices())
		{
			ExtendTestCode extendCode(pC.kernelParameters, device);
			fill_buffers();
			extendCode.runExtendKernel(d_out, d_in, d_rand, NUM_ELEMENTS);
			verify();
		}
		BOOST_TEST_MESSAGE("Tested CPU");
	}
	void verify();
	void verify(hmc_complex, hmc_complex);
	void fill_buffers();
	Matrixsu2 * h_in;
	Matrixsu3 * h_out;
	cl_int * h_rand;
	hardware::buffers::Plain<Matrixsu2> * d_in;
	hardware::buffers::Plain<Matrixsu3> * d_out;
	hardware::buffers::Plain<cl_int> * d_rand;
};

BOOST_AUTO_TEST_CASE( EXTEND )
{
	GaugefieldTestParameters parametersForThisTest {LatticeExtents{4,4}, GaugefieldFillType::cold};
	hardware::HardwareParametersMockup hardwareParameters(parametersForThisTest.ns,parametersForThisTest.nt);
	hardware::code::OpenClKernelParametersMockup kernelParameters(parametersForThisTest.ns,parametersForThisTest.nt);
	ParameterCollection parameterCollection(hardwareParameters, kernelParameters);
	ExtendTester tester(parameterCollection, parametersForThisTest);
}

void ExtendTester::verify(hmc_complex left, hmc_complex right)
{
	BOOST_REQUIRE_EQUAL(left.re, right.re);
	BOOST_REQUIRE_EQUAL(left.im, right.im);
}

void ExtendTester::verify()
{
	// get stuff from device
	d_out->dump(h_out);

	for(size_t i = 0; i < NUM_ELEMENTS; ++i) {
		const int rand = h_rand[i];
		const Matrixsu3 m = h_out[i];

		verify(m.e00, (rand == 2) ? hmc_complex_one : hmc_complex_zero);
		verify(m.e01, hmc_complex_zero);
		verify(m.e02, hmc_complex_zero);
		verify(m.e10, hmc_complex_zero);
		verify(m.e11, (rand == 3) ? hmc_complex_one : hmc_complex_zero);
		verify(m.e12, hmc_complex_zero);
		verify(m.e20, hmc_complex_zero);
		verify(m.e21, hmc_complex_zero);
		verify(m.e22, (rand == 1) ? hmc_complex_one : hmc_complex_zero);
	}
}

void ExtendTester::fill_buffers()
{
	const hardware::Device * device = code->get_device();

	h_in = new Matrixsu2[NUM_ELEMENTS];
	BOOST_REQUIRE(h_in);
	for(int i = 0; i < NUM_ELEMENTS; ++i) {
		h_in[i].e00 = hmc_complex_zero;
		h_in[i].e01 = hmc_complex_zero;
		h_in[i].e10 = hmc_complex_zero;
		h_in[i].e11 = hmc_complex_zero;
	}

	d_in = new hardware::buffers::Plain<Matrixsu2>(NUM_ELEMENTS, device);
	d_in->load(h_in);

	// for simplicity initialize input with 0

	h_out = new Matrixsu3[NUM_ELEMENTS];
	BOOST_REQUIRE(h_out);
	d_out = new hardware::buffers::Plain<Matrixsu3>(NUM_ELEMENTS, device);

	h_rand = new cl_int[NUM_ELEMENTS];
	BOOST_REQUIRE(h_rand);
	for(int i = 0; i < NUM_ELEMENTS; ++i) {
		h_rand[i] = (i % 3) + 1; // high quality random numbers ;)
	}
	d_rand = new hardware::buffers::Plain<cl_int>(NUM_ELEMENTS, device);
	d_rand->load(h_rand);
}

