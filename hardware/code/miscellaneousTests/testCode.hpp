/*
 * Copyright 2015 Christopher Pinke
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

# pragma once

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
		testKernel = 0;
		err = 0;
	};

	virtual ~TestCode()
	{
		if(testKernel)
			clReleaseKernel(testKernel);
	};

	virtual void runTestKernel(const hardware::buffers::SU3 *, const hardware::buffers::Plain<hmc_float> *, const int, const int) {};

	cl_kernel testKernel;
	cl_int err;
};

struct OtherKernelTester : public GaugefieldTester
{
	OtherKernelTester(const std::string kernelName, const ParameterCollection pC, const GaugefieldTestParameters tP, const ReferenceValues rV):
		GaugefieldTester(kernelName, pC, tP, rV), testCode(0), tP(tP)
	{
		GaugefieldCreator gf(tP.latticeExtents);
		gaugefieldBuffer = new hardware::buffers::SU3(calculateGaugefieldSize(tP.latticeExtents), this->device);
		const Matrixsu3 * gf_host = gf.createGaugefield(tP.fillType);
        device->getGaugefieldCode()->importGaugefield(gaugefieldBuffer, gf_host);
        delete[] gf_host;
		out = new hardware::buffers::Plain<hmc_float> (LatticeExtents(tP.latticeExtents).getLatticeVolume(), device);

		if(device->get_device_type() == CL_DEVICE_TYPE_GPU)
		{
			gs = LatticeExtents(tP.latticeExtents).getLatticeVolume();
			ls = 64;
		}
		else
		{
			gs = device->get_num_compute_units();
			ls = 1;
		}
	}
	~OtherKernelTester()
	{
		hmc_float * host_out = new hmc_float[LatticeExtents(tP.latticeExtents).getLatticeVolume()];
		out->dump(host_out);

		hmc_float result = 0;
		for(unsigned int i = 0; i < LatticeExtents(tP.latticeExtents).getLatticeVolume(); i++) {
			result += host_out[i];
		}
		kernelResult.at(0) = result;

		if(testCode)
			delete testCode;
		if(host_out)
			delete host_out;
		if(out)
			delete out;
	}
	TestCode * testCode;
	int ls, gs;
	const GaugefieldTestParameters tP;
	hardware::buffers::Plain<hmc_float> * out;
//protected:
	const hardware::buffers::SU3 * gaugefieldBuffer;
};



