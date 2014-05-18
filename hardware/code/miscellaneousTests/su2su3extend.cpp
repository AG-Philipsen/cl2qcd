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

#include "../hardware/system.hpp"
#include "../hardware/device.hpp"
#include "../hardware/code/opencl_module.hpp"
#include "../hardware/code/gaugefield.hpp"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE SU2 SU3 Extend
#include <boost/test/unit_test.hpp>

#define NUM_ELEMENTS 1024
#define LOCAL_SIZE 128

class Device : public hardware::code::Opencl_Module {

	cl_kernel extendKernel;

	void fill_kernels();
	void clear_kernels();
protected:
	virtual size_t get_read_write_size(const std::string&) const {
		return 0;
	};
	virtual uint64_t get_flop_size(const std::string&) const {
		return 0;
	};
public:
	Device(const meta::Inputparameters& params, hardware::Device * device) : Opencl_Module(params, device) {
		fill_kernels();
	};
	~Device() {
		clear_kernels();
	};

	void runExtendKernel(const hardware::buffers::Plain<Matrixsu3> * out, const hardware::buffers::Plain<Matrixsu2> * in, const hardware::buffers::Plain<cl_int> * d_rand, cl_ulong elems);
};

class Test {

public:
	Test(const meta::Inputparameters& params, hardware::Device* device)
		: params(params), code(params, device) {
		fill_buffers();
	};
	~Test() {
		clear_buffers();
	}

	void verify();
	void runExtendKernel();

private:
	void verify(hmc_complex, hmc_complex);
	void fill_buffers();
	void clear_buffers();
	Matrixsu2 * h_in;
	Matrixsu3 * h_out;
	cl_int * h_rand;
	hardware::buffers::Plain<Matrixsu2> * d_in;
	hardware::buffers::Plain<Matrixsu3> * d_out;
	hardware::buffers::Plain<cl_int> * d_rand;

	const meta::Inputparameters& params;
	Device code;
};

BOOST_AUTO_TEST_CASE( CPU )
{
	const char* _params_cpu[] = {"foo", "--use_gpu=false"};
	meta::Inputparameters params_cpu(2, _params_cpu);
	hardware::System system(params_cpu);
for(auto device: system.get_devices()) {
		Test dummy(params_cpu, device);
		dummy.runExtendKernel();
		dummy.verify();
	}
	BOOST_MESSAGE("Tested CPU");
}

BOOST_AUTO_TEST_CASE( GPU )
{
	const char* _params_gpu[] = {"foo", "--use_gpu=true"};
	meta::Inputparameters params_gpu(2, _params_gpu);
	hardware::System system(params_gpu);
for(auto device: system.get_devices()) {
		Test dummy(params_gpu, device);
		dummy.runExtendKernel();
		dummy.verify();
	}
	BOOST_MESSAGE("Tested GPU");
}

void Test::fill_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	hardware::Device * device = code.get_device();

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

void Device::fill_kernels()
{
	extendKernel = createKernel("extendKernel") << get_device()->get_gaugefield_code()->get_sources() << "../hardware/code/miscellaneousTests/su2su3extend.cl";
}

void Test::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	delete d_in;
	delete d_out;
	delete d_rand;

	delete[] h_in;
	delete[] h_out;
	delete[] h_rand;
}

void Device::clear_kernels()
{
	clReleaseKernel(extendKernel);
}

void Device::runExtendKernel(const hardware::buffers::Plain<Matrixsu3> * out, const hardware::buffers::Plain<Matrixsu2> * in, const hardware::buffers::Plain<cl_int> * d_rand, cl_ulong elems)
{
	cl_int err;
	err = clSetKernelArg(extendKernel, 0, sizeof(cl_mem), out->get_cl_buffer());
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(extendKernel, 1, sizeof(cl_mem), in->get_cl_buffer());
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(extendKernel, 2, sizeof(cl_mem), d_rand->get_cl_buffer());
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(extendKernel, 3, sizeof(cl_ulong), &elems);
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	get_device()->enqueue_kernel(extendKernel, NUM_ELEMENTS);
}

void Test::verify(hmc_complex left, hmc_complex right)
{
	BOOST_REQUIRE_EQUAL(left.re, right.re);
	BOOST_REQUIRE_EQUAL(left.im, right.im);
}


void Test::verify()
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

void Test::runExtendKernel()
{
	code.runExtendKernel(d_out, d_in, d_rand, NUM_ELEMENTS);
}
