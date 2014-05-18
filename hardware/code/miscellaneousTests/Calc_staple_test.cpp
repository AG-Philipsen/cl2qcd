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

#include "../meta/util.hpp"
#include "../physics/lattices/gaugefield.hpp"
#include "../hardware/device.hpp"

#include "../../../physics/observables/gaugeObservables.h"

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE staple_test
#include <boost/test/unit_test.hpp>

std::string const exec_name = "staple_test";

class Code : public hardware::code::Opencl_Module {

	cl_kernel testKernel;
protected:
	virtual size_t get_read_write_size(const std::string&) const {
		return 0;
	};
	virtual uint64_t get_flop_size(const std::string&) const {
		return 0;
	};
public:
	Code(const meta::Inputparameters& params, hardware::Device * device) : Opencl_Module(params, device) {
	  logger.info() << "init code";
	  fill_kernels();
	  logger.info() << "code done";
	};
	~Code() {
		clear_kernels();
	};

	void runTestKernel(const hardware::buffers::SU3 * gf, const hardware::buffers::Plain<hmc_float> * out, int gs, int ls);
	void fill_kernels();
	void clear_kernels();
};


class Dummyfield {
public:
	Dummyfield(const hardware::System& system) : device(system.get_devices().at(0)), params(system.get_inputparameters()), code(params, device), prng(system), gf(system, prng) {
		meta::print_info_hmc(system.get_inputparameters());
		fill_buffers();
	};
	~Dummyfield() {
		clear_buffers();
	}
	hmc_float get_squarenorm();
	hmc_float runTestKernel();

private:
	void fill_buffers();
	void clear_buffers();
	const hardware::buffers::Plain<hmc_float> * out;
	hmc_float * host_out;
	hardware::Device * const device;
	const meta::Inputparameters& params;
	Code code;
	physics::PRNG prng;
public:
	const physics::lattices::Gaugefield gf;
};

void Dummyfield::fill_buffers()
{
	// don't invoke parent function as we don't require the original buffers
	int NUM_ELEMENTS = meta::get_vol4d(params);
	host_out = new hmc_float[NUM_ELEMENTS];
	BOOST_REQUIRE(host_out);

	out = new hardware::buffers::Plain<hmc_float>(NUM_ELEMENTS, device);
}

void Code::fill_kernels()
{
  logger.info() << "create Kernel";
	testKernel = createKernel("staple_test") << get_device()->get_gaugefield_code()->get_sources() << "../hardware/code/miscellaneousTests/staple_test.cl";
	logger.info() << "done";
}

void Dummyfield::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers
	delete out;
	delete[] host_out;
}

void Code::clear_kernels()
{
	clReleaseKernel(testKernel);
}

void Code::runTestKernel(const hardware::buffers::SU3 * gf, const hardware::buffers::Plain<hmc_float> * out, int gs, int ls)
{
  logger.info() << "call kernel";
	cl_int err;
	err = clSetKernelArg(testKernel, 0, sizeof(cl_mem), gf->get_cl_buffer());
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);
	err = clSetKernelArg(testKernel, 1, sizeof(cl_mem), out->get_cl_buffer());
	BOOST_REQUIRE_EQUAL(CL_SUCCESS, err);

	get_device()->enqueue_kernel(testKernel, gs, ls);
}

hmc_float Dummyfield::runTestKernel()
{
	hmc_float res = 0;
	int gs, ls;
	if(device->get_device_type() == CL_DEVICE_TYPE_GPU) {
		gs = meta::get_vol4d(params);
		ls = 64;
	} else {
		gs = device->get_num_compute_units();
		ls = 1;
	}
	code.runTestKernel(gf.get_buffers()[0], out, gs, ls);

	int NUM_ELEMENTS = meta::get_vol4d(params);
	//copy the result of the kernel to host
	out->dump(host_out);

	//sum up all elements in the result buffer
	for(int i = 0; i < NUM_ELEMENTS; i++) {
		res += host_out[i];
	}
	return res;
}


BOOST_AUTO_TEST_CASE( STAPLE_TEST )
{
	logger.info() << "Test kernel";
	logger.info() << "\tcalc_staple";
	logger.info() << "against reference value";

	int param_expect = 4;
	logger.info() << "expect parameters:";
	logger.info() << "\texec_name\tinputfile\tgpu_usage\trec12_usage";
	//get number of parameters
	int num_par = boost::unit_test::framework::master_test_suite().argc;
	if(num_par < param_expect) {
		logger.fatal() << "need more inputparameters! Got only " << num_par << ", expected " << param_expect << "! Aborting...";
		exit(-1);
	}

	//get input file that has been passed as an argument
	const char*  inputfile =  boost::unit_test::framework::master_test_suite().argv[1];
	logger.info() << "inputfile used: " << inputfile;
	//get use_gpu = true/false that has been passed as an argument
	const char*  gpu_opt =  boost::unit_test::framework::master_test_suite().argv[2];
	logger.info() << "GPU usage: " << gpu_opt;
	//get use_rec12 = true/false that has been passed as an argument
	const char* rec12_opt =  boost::unit_test::framework::master_test_suite().argv[3];
	logger.info() << "rec12 usage: " << rec12_opt;

	logger.info() << "Init device";
	const char* _params_cpu[] = {"foo", inputfile, gpu_opt, rec12_opt, "--device=0"};
	meta::Inputparameters params(param_expect + 1, _params_cpu);
	hardware::System system(params);
	Dummyfield cpu(system);
	logger.info() << "gaugeobservables: ";
	physics::gaugeObservables obs(&params );
	obs.measureGaugeObservables(&cpu.gf, 0);
	logger.info() << "Run kernel";
	logger.info() << "running test kernel";
	hmc_float cpu_res = cpu.runTestKernel();
	logger.info() << "result:";
	logger.info() << cpu_res;

	logger.info() << "Choosing reference value and acceptance precision";
	hmc_float ref_val = params.get_test_ref_value();
	logger.info() << "reference value:\t" << ref_val;
	hmc_float prec = params.get_solver_prec();
	logger.info() << "acceptance precision: " << prec;

	logger.info() << "Compare result to reference value";
	BOOST_REQUIRE_CLOSE(cpu_res, ref_val, prec);
	logger.info() << "Done";
	BOOST_MESSAGE("Test done");
}

