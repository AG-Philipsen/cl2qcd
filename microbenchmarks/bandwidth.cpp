/** @file
 * A simple microbenchmark to check device bandwidth for different
 * numbers of threada and streaming.
 *
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

#include <string>
#include <stdexcept>

#include <boost/program_options.hpp>

#include "../hardware/system.hpp"
#include "../hardware/device.hpp"
#include "../hardware/code/gaugefield.hpp"
#include "../host_functionality/logger.hpp"
#include "../executables/exceptions.h"
#include "../klepsydra/klepsydra.hpp"

#include "../common_header_files/types_fermions.h"


namespace po = boost::program_options;

/**
 * Selector type for the base type of the copy operations.
 */
enum copyType {
  type_invalid,
  type_float,
  type_su3,
  type_su3SOA,
  type_su3SOcplxA,
  type_spinor,
  type_spinorSOA,
  type_spinorSOcplxA,
  type_spinorLocal,
  type_spinorSOApy
};

size_t getTypeSize(copyType type);

class Device : public hardware::code::Opencl_Module {

private:
	cl_kernel floatKernel;
	cl_kernel su3Kernel;
	cl_kernel su3SOAKernel;
	cl_kernel su3SOcplxAKernel;
	cl_kernel spinorKernel;
	cl_kernel spinorSOAKernel;
	cl_kernel spinorSOcplxAKernel;
	cl_kernel spinorLocalKernel;
	cl_kernel spinorSOApyKernel;

	template<typename T> void runKernel(size_t groups, cl_ulong threads_per_group, cl_ulong elems, cl_kernel kernel, const hardware::buffers::Plain<cl_char> * in, const hardware::buffers::Plain<cl_char> * out);

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
	virtual ~Device() {
		clear_kernels();
	};

	void runKernel(copyType copy_type, size_t groups, cl_ulong threads_per_group, cl_ulong elems, const hardware::buffers::Plain<cl_char> * in, const hardware::buffers::Plain<cl_char> * out);
};

class Test {

public:
	Test(const hardware::System& system, size_t maxMemSize)
		: maxMemSize(maxMemSize), device(system.get_devices().at(0)), code(system.get_inputparameters(), device) {
		fill_buffers();
	};

	void runKernel(copyType copy_type, size_t groups, cl_ulong threads_per_group, cl_ulong elems);

	~Test();
private:
	void verify(hmc_complex, hmc_complex);
	void fill_buffers();
	const hardware::buffers::Plain<cl_char> * in, * out;
	size_t maxMemSize;
	hardware::Device * const device;
	Device code;
};

int main(int argc, char** argv)
{
	po::options_description desc("Allowed options");
	desc.add_options()
	("help,h", "Produce this help message")
	("elements,e", po::value<cl_ulong>()->default_value(100000), "How many elements to use.") // conflicts with single
	("stepelements,se", po::value<cl_ulong>(), "Step size for element count sweeping")
	("threads,t", po::value<cl_ulong>()->default_value(64), "The number of threads to use per groups")
	("groups,g", po::value<cl_ulong>(), "Specify a fixed number of groups.")
	("stepthreads,st", po::value<cl_ulong>()->default_value(0), "Step size for thread per group sweeping")
	("single", "Copy only a single element per thread")
	("type,d", po::value<std::string>()->default_value("float"), "The data type to copy");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	// Check if help is needed
	if( vm.count( "help" ) ) {
		std::cout << desc << '\n';
		return 0;
	}

	// parse type
	std::map<std::string, copyType> type_map;
	type_map["float"] = type_float;
	type_map["su3"] = type_su3;
	type_map["su3SOA"] = type_su3SOA;
	type_map["su3SOcplxA"] = type_su3SOcplxA;
	type_map["spinor"] = type_spinor;
	type_map["spinorSOA"] = type_spinorSOA;
	type_map["spinorSOApy"] = type_spinorSOApy;
	type_map["spinorSOcplxA"] = type_spinorSOcplxA;
	type_map["spinorLocal"] = type_spinorLocal;

	const copyType copy_type = type_map[vm["type"].as<std::string>()];
	if(!copy_type) {
		logger.error() << "Please select one of the following types: float(default), su3, su3SOA, su3SOcplxA, spinor, spinorSOA, spinorSOApy, spinorSOcplxA, spinorLocal";
		return 1;
	}
	logger.info() << "Using " << vm["type"].as<std::string>() << " as load/store datatype";

	if(vm.count("single") && vm.count("elements")) {
		logger.error() << "You can either use one element per thread or define the global number of elements";
	}

	meta::Inputparameters params(0, 0);
	hardware::System system(params, true);

	if(vm.count("stepelements")) {
		logger.info() << "Sweeping element count for fixed thread count.";
		cl_ulong max_elements = vm["elements"].as<cl_ulong>();
		cl_ulong step_elements = vm["stepelements"].as<cl_ulong>();
		cl_ulong threads = vm["threads"].as<cl_ulong>();
		cl_ulong groups;
		if(vm.count("groups")) {
			groups = vm["groups"].as<cl_ulong>();
		} else {
			groups = 20 * 8; // taken out of the air and tuned for Cypress
		}
		if(vm.count("single")) {
			logger.fatal() << "Single element per thread mode has not been implemented in element count sweeping mod";
		} else {
			Test dev(system, max_elements * getTypeSize(copy_type));
			size_t elements = 1;
			dev.runKernel(copy_type, groups, threads, elements);
			for(elements = step_elements; elements <= max_elements; elements += step_elements) {
				dev.runKernel(copy_type, groups, threads, elements);
			}
		}

	} else if(vm.count("groups")) {
		cl_ulong max_threads = vm["threads"].as<cl_ulong>();
		cl_ulong groups = vm["groups"].as<cl_ulong>();
		cl_ulong step_threads = vm["stepthreads"].as<cl_ulong>();
		cl_ulong min_threads;
		if(step_threads) {
			logger.info() << "Scanning number active threads per group required for maximum memory throughput";
			min_threads = step_threads;
		} else {
			min_threads = max_threads;
			step_threads = max_threads + 1;
		}
		if(vm.count("single")) {
			logger.info() << "Using a single element per thread";
			Test dev(system, groups * max_threads * getTypeSize(copy_type));
			if(step_threads <= max_threads) {
				size_t threads = 1;
				dev.runKernel(copy_type, groups, threads, groups * threads);
			}
			for(size_t threads = min_threads; threads <= max_threads; threads += step_threads) {
				dev.runKernel(copy_type, groups, threads, groups * threads);
			}
		} else {
			const cl_ulong elems = vm["elements"].as<cl_ulong>();
			logger.info() << "Keeping number of elements fixed at " << elems;
			Test dev(system, elems * getTypeSize(copy_type));
			if(step_threads <= max_threads) {
				size_t threads = 1;
				dev.runKernel(copy_type, groups, threads, groups * threads);
			}
			for(size_t threads = min_threads; threads <= max_threads; threads += step_threads) {
				dev.runKernel(copy_type, groups, threads, elems);
			}
		}
	} else {
		logger.info() << "Scanning number of wavefronts required for maximum memory throughput";
		cl_ulong threads = vm["threads"].as<cl_ulong>();
		cl_ulong max_groups = 20 * 8; // the value comes out of nowhere and is tuned for Cypress

		if(vm.count("single")) {
			logger.info() << "Using a single element per thread";
			Test dev(system, max_groups * threads * getTypeSize(copy_type));
			for(size_t groups = 1; groups <= max_groups; ++groups) {
				dev.runKernel(copy_type, groups, threads, groups * threads);
			}
		} else {
			const cl_ulong elems = vm["elements"].as<cl_ulong>();
			logger.info() << "Keeping number of elements fixed at " << elems;
			Test dev(system, elems * getTypeSize(copy_type));
			for(size_t groups = 1; groups <= max_groups; ++groups) {
				dev.runKernel(copy_type, groups, threads, elems);
			}
		}
	}

	logger.info() << "Done";

	return 0;
}

void Test::fill_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	// kernels might stride, always overallocate memory
	// e.g. spinor might pad up to 24 * 8 KiB = 128 KiB -> 1 MiB should be save
	size_t pad_buf = 1024 * 1024;
	size_t allocMemSize = maxMemSize + pad_buf;

	logger.info() << "Allocating buffers of " << allocMemSize << " bytes.";

	in = new hardware::buffers::Plain<cl_char>(allocMemSize, device);
	out = new hardware::buffers::Plain<cl_char>(allocMemSize, device);
}

void Device::fill_kernels()
{
	ClSourcePackage basic_opencl_code = get_device()->getGaugefieldCode()->get_sources();
	floatKernel = createKernel("copyFloat") << basic_opencl_code << "types_fermions.h" << "../microbenchmarks/bandwidth.cl";
	su3Kernel = createKernel("copySU3") << basic_opencl_code << "types_fermions.h" << "../microbenchmarks/bandwidth.cl";
	su3SOAKernel = createKernel("copySU3SOA") << basic_opencl_code << "types_fermions.h" << "../microbenchmarks/bandwidth.cl";
	su3SOcplxAKernel = createKernel("copySU3SOcplxA") << basic_opencl_code << "types_fermions.h" << "../microbenchmarks/bandwidth.cl";
	spinorKernel = createKernel("copySpinor") << basic_opencl_code << "types_fermions.h" << "../microbenchmarks/bandwidth.cl";
	spinorSOAKernel = createKernel("copySpinorSOA") << basic_opencl_code << "types_fermions.h" << "../microbenchmarks/bandwidth.cl";
	spinorSOApyKernel = createKernel("copyDpSpinorFullestSOARestricted") << basic_opencl_code << "types_fermions.h" << "../microbenchmarks/bandwidth_spinorSOApy.cl";
	spinorSOcplxAKernel = createKernel("copySpinorSOcplxA") << basic_opencl_code << "types_fermions.h" << "../microbenchmarks/bandwidth.cl";
	spinorLocalKernel = createKernel("copySpinorLocal") << basic_opencl_code << "types_fermions.h" << "../microbenchmarks/bandwidth.cl";
}

Test::~Test()
{
	// don't invoke parent function as we don't require the original buffers
	delete in;
	delete out;
}

void Device::clear_kernels()
{
	clReleaseKernel(floatKernel);
	clReleaseKernel(su3Kernel);
	clReleaseKernel(su3SOAKernel);
	clReleaseKernel(su3SOcplxAKernel);
	clReleaseKernel(spinorKernel);
	clReleaseKernel(spinorSOAKernel);
	clReleaseKernel(spinorSOApyKernel);
	clReleaseKernel(spinorSOcplxAKernel);
	clReleaseKernel(spinorLocalKernel);
}


void Device::runKernel(copyType copy_type, size_t groups, cl_ulong threads_per_group, cl_ulong elems, const hardware::buffers::Plain<cl_char> * in, const hardware::buffers::Plain<cl_char> * out)
{
	switch(copy_type) {
		case type_float:
			runKernel<hmc_float>(groups, threads_per_group, elems, floatKernel, in, out);
			return;
		case type_su3:
			runKernel<Matrixsu3>(groups, threads_per_group, elems, su3Kernel, in, out);
			return;
		case type_su3SOA:
			runKernel<Matrixsu3>(groups, threads_per_group, elems, su3SOAKernel, in, out);
			return;
		case type_su3SOcplxA:
			runKernel<Matrixsu3>(groups, threads_per_group, elems, su3SOcplxAKernel, in, out);
			return;
		case type_spinor:
			runKernel<spinor>(groups, threads_per_group, elems, spinorKernel, in, out);
			return;
		case type_spinorSOA:
			runKernel<spinor>(groups, threads_per_group, elems, spinorSOAKernel, in, out);
			return;
		case type_spinorSOApy:
			runKernel<spinor>(groups, threads_per_group, elems, spinorSOApyKernel, in, out);
			return;
		case type_spinorSOcplxA:
			runKernel<spinor>(groups, threads_per_group, elems, spinorSOcplxAKernel, in, out);
			return;
		case type_spinorLocal:
			runKernel<spinor>(groups, threads_per_group, elems, spinorLocalKernel, in, out);
			return;
		default:
			throw std::invalid_argument("runKernel has not been implemented for this type");
	}
}

template<typename T> void Device::runKernel(size_t groups, cl_ulong threads_per_group, cl_ulong elems, cl_kernel kernel, const hardware::buffers::Plain<cl_char> * in, const hardware::buffers::Plain<cl_char> * out)
{
	cl_int err = CL_SUCCESS;

	// TODO adjust threads_per_group for kernel invocation to proper size but keep requested value for kernel arg
	size_t local_threads = threads_per_group;
	size_t total_threads = groups * local_threads;

	err = clSetKernelArg(kernel, 0, sizeof(cl_mem), out->get_cl_buffer());
	if(err) {
		logger.fatal() << "Failed to set kernel argument: " << err;
		throw Opencl_Error(err);
	}
	err = clSetKernelArg(kernel, 1, sizeof(cl_mem), in->get_cl_buffer());
	if(err) {
		logger.fatal() << "Failed to set kernel argument: " << err;
		throw Opencl_Error(err);
	}
	err = clSetKernelArg(kernel, 2, sizeof(cl_ulong), &elems);
	if(err) {
		logger.fatal() << "Failed to set kernel argument: " << err;
		throw Opencl_Error(err);
	}
	err = clSetKernelArg(kernel, 3, sizeof(cl_ulong), &threads_per_group);
	if(err) {
		logger.fatal() << "Failed to set kernel argument: " << err;
		throw Opencl_Error(err);
	}

	size_t num_meas = 10;

	hardware::Device * dev = get_device();
	dev->enqueue_kernel(kernel, total_threads, local_threads);
	dev->synchronize();
	klepsydra::Monotonic timer;
	for(size_t i = 0; i < num_meas; ++i)
		dev->enqueue_kernel(kernel, total_threads, local_threads);
	dev->synchronize();
	int64_t kernelTime = timer.getTime() / num_meas;
	if(err) {
		logger.fatal() << "Failed to execute kernel: " << err;
		throw Opencl_Error(err);
	}

	cl_ulong eventTime = dev->getProfilingData(kernel).get_total_time() / num_meas;

	// format is: #groups #threads per group #elements #copied memory in bytes #copy time in mus #bandwidth in megabytes
	// FIXME sizeof can give broken results in case of aligned types (gross size not equal to net content size)
	std::cout << groups * threads_per_group << ' ' << groups << ' ' << threads_per_group << ' ' << elems << ' ' << elems * sizeof(T) << ' ' << kernelTime << ' ' << (2 * elems * sizeof(T) / eventTime)   << ' ' << eventTime << std::endl;
}

void Test::runKernel(copyType copy_type, size_t groups, cl_ulong threads_per_group, cl_ulong elems)
{
	code.runKernel(copy_type, groups, threads_per_group, elems, in, out);
}

size_t getTypeSize(copyType type)
{
	switch(type) {
		case type_float:
			return sizeof(hmc_float);
		case type_su3:
		case type_su3SOA:
		case type_su3SOcplxA:
			return sizeof(Matrixsu3);
		case type_spinor:
		case type_spinorSOcplxA:
		case type_spinorSOA:
		case type_spinorSOApy:
		case type_spinorLocal:
			return sizeof(spinor);
		default:
			throw std::invalid_argument("getTypeSize has not been implemented for this type.");
	}
}
