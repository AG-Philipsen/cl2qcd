/** @file
 * A simple microbenchmark to check device bandwidth for different
 * numbers of threada and streaming.
 */

#include <string>

#include <boost/program_options.hpp>

#include "../host_random.h"
#include "../opencl.h"
#include "../logger.hpp"

namespace po = boost::program_options;

extern std::string const version;
std::string const version = "0.1";

Random rnd(15);

const size_t MAX_MEM_SIZE = 16 * 1024 * 1024;

class Device : public Opencl {

private:
	inputparameters params;
	usetimer timer;
	cl_kernel floatKernel;
	cl_kernel su3Kernel;
	cl_mem in, out;

	template<typename T> void runKernel(size_t groups, cl_ulong threads_per_group, cl_ulong elems, cl_kernel kernel);

public:
	Device(cl_device_type device_type) : Opencl() {
		Opencl::init(device_type, &timer, &params); /* init in body for proper this-pointer */
	};
	virtual hmc_error fill_kernels_file();
	virtual hmc_error fill_buffers();
	virtual hmc_error fill_kernels();
	virtual hmc_error clear_buffers();
	virtual hmc_error clear_kernels();
	~Device() {
		finalize();
	};

	void runFloatKernel(size_t groups, cl_ulong threads_per_group, cl_ulong elems);
	void runSU3Kernel(size_t groups, cl_ulong threads_per_group, cl_ulong elems);
};

int main(int argc, char** argv)
{
	po::options_description desc("Allowed options");
	desc.add_options()
	("help,h", "Produce this help message")
//	     ("elements,e", po::value<cl_ulong>()->default_value(100000), "How many elements to use.") // conflicts with single
	("threads,t", po::value<cl_ulong>()->default_value(64), "The number of threads to use per groups (maximum if groups is set)")
	("groups,g", po::value<cl_ulong>(), "Vary number of threads per group for a fixed number of groups") // default is to vary number of groups for 64 threads per groups
	("single", "Copy only a single element per thread")
	("su3", "Use SU3 datastructure instead of float");

	po::variables_map vm;
	po::store(po::parse_command_line(argc, argv, desc), vm);
	po::notify(vm);

	// Check if help is needed
	if( vm.count( "help" ) ) {
		std::cout << desc << '\n';
		return 0;
	}

	Device dev(CL_DEVICE_TYPE_GPU);

	const bool use_su3 = vm.count("su3");
	if(use_su3)
		logger.info() << "Using su3 matrices as load/store datatype";
	else
		logger.info() << "Using floating point data type";

	if(vm.count("groups")) {
		logger.info() << "Scanning number active threads per group required for maximum memory throughput";
		cl_ulong max_threads = vm["threads"].as<cl_ulong>();
		cl_ulong groups = vm["groups"].as<cl_ulong>();

		if(vm.count("single")) {
			logger.info() << "Using a single element per thread";
			for(size_t threads = 1; threads <= max_threads; ++threads) {
				if(use_su3)
					dev.runSU3Kernel(groups, threads, groups * threads);
				else
					dev.runFloatKernel(groups, threads, groups * threads);
			}
		} else {
			const cl_ulong elems = MAX_MEM_SIZE / (use_su3 ? sizeof(Matrixsu3) : sizeof(hmc_float) );
			logger.info() << "Keeping number of elements fixed at " << elems;
			for(size_t threads = 1; threads <= max_threads; ++threads) {
				if(use_su3)
					dev.runSU3Kernel(groups, threads, elems);
				else
					dev.runFloatKernel(groups, threads, elems);
			}
		}

	} else {
		logger.info() << "Scanning number of wavefronts required for maximum memory throughput";
		cl_ulong threads = vm["threads"].as<cl_ulong>();
		cl_ulong max_groups = 20 * 8; // the value comes out of nowhere and is tuned for Cypress

		if(vm.count("single")) {
			logger.info() << "Using a single element per thread";
			for(size_t groups = 1; groups <= max_groups; ++groups) {
				if(use_su3)
					dev.runSU3Kernel(groups, threads, groups * threads);
				else
					dev.runFloatKernel(groups, threads, groups * threads);
			}
		} else {
			const cl_ulong elems = MAX_MEM_SIZE / (use_su3 ? sizeof(Matrixsu3) : sizeof(hmc_float) );
			logger.info() << "Keeping number of elements fixed at " << elems;
			for(size_t groups = 1; groups <= max_groups; ++groups) {
				if(use_su3)
					dev.runSU3Kernel(groups, threads, elems);
				else
					dev.runFloatKernel(groups, threads, elems);
			}
		}
	}

	logger.info() << "Done";

	return 0;
}

hmc_error Device::fill_kernels_file()
{
	//give a list of all kernel-files
	Opencl::fill_kernels_file();

	cl_kernels_file.push_back("microbenchmarks/bandwidth.cl");

	logger.debug() << "Filled kernels file list";

	return HMC_SUCCESS;
}

hmc_error Device::fill_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	cl_int err;

	in = clCreateBuffer(context, CL_MEM_READ_ONLY, MAX_MEM_SIZE, 0, &err );
	if(err) {
		logger.fatal() << "Unable to allocate memory on device";
		exit(HMC_OCLERROR);
	}

	out = clCreateBuffer(context, CL_MEM_WRITE_ONLY, MAX_MEM_SIZE, 0, &err );
	if(err) {
		logger.fatal() << "Unable to allocate memory on device";
		exit(HMC_OCLERROR);
	}

	return HMC_SUCCESS;
}

hmc_error Device::fill_kernels()
{
	int clerr = HMC_SUCCESS;

	// don't invoke parent function as we don't require the original kernels

	floatKernel = clCreateKernel(clprogram, "copyFloat", &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "...creating float kernel failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	if( logger.beDebug() )
		printResourceRequirements( floatKernel );

	su3Kernel = clCreateKernel(clprogram, "copySU3", &clerr);
	if(clerr != CL_SUCCESS) {
		cout << "...creating su3 kernel failed, aborting." << endl;
		exit(HMC_OCLERROR);
	}
	if( logger.beDebug() )
		printResourceRequirements( su3Kernel );
	return clerr;
}

hmc_error Device::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	clReleaseMemObject(in);
	clReleaseMemObject(out);

	return HMC_SUCCESS;
}

hmc_error Device::clear_kernels()
{
	int clerr = HMC_SUCCESS;

	// don't invoke parent function as we don't require the original kernels

	clReleaseKernel(floatKernel);
	clReleaseKernel(su3Kernel);

	return clerr;
}


void Device::runFloatKernel(size_t groups, cl_ulong threads_per_group, cl_ulong elems)
{
	runKernel<hmc_float>(groups, threads_per_group, elems, floatKernel);
}

void Device::runSU3Kernel(size_t groups, cl_ulong threads_per_group, cl_ulong elems)
{
	runKernel<Matrixsu3>(groups, threads_per_group, elems, su3Kernel);
}

template<typename T> void Device::runKernel(size_t groups, cl_ulong threads_per_group, cl_ulong elems, cl_kernel kernel)
{
	cl_int err = CL_SUCCESS;

	// make sure memory buffers are large enough
	if( elems * sizeof(T) > MAX_MEM_SIZE ) {
		logger.fatal() << "Not enough memory for the given number of elements";
		exit(-1);
	}

	// TODO adjust threads_per_group for kernel invocation to proper size but keep requested value for kernel arg
	size_t local_threads = threads_per_group;
	size_t total_threads = groups * local_threads;

	err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &in);
	if(err) {
		logger.fatal() << "Failed to set kernel argument: " << err;
		exit(HMC_OCLERROR);
	}
	err = clSetKernelArg(kernel, 1, sizeof(cl_mem), &in);
	if(err) {
		logger.fatal() << "Failed to set kernel argument: " << err;
		exit(HMC_OCLERROR);
	}
	err = clSetKernelArg(kernel, 2, sizeof(cl_ulong), &elems);
	if(err) {
		logger.fatal() << "Failed to set kernel argument: " << err;
		exit(HMC_OCLERROR);
	}
	err = clSetKernelArg(kernel, 3, sizeof(cl_ulong), &threads_per_group);
	if(err) {
		logger.fatal() << "Failed to set kernel argument: " << err;
		exit(HMC_OCLERROR);
	}

	size_t num_meas = 10;

	klepsydra::Monotonic timer;

	enqueueKernel(kernel, total_threads, local_threads);
	err = clFinish(queue);
	if(err) {
		logger.fatal() << "Failed to execute kernel: " << err;
		exit(HMC_OCLERROR);
	}
	for(size_t i = 1; i < num_meas; ++i)
		enqueueKernel(kernel, total_threads, local_threads);
	err = clFinish(queue);
	if(err) {
		logger.fatal() << "Failed to execute kernel: " << err;
		exit(HMC_OCLERROR);
	}
	int64_t kernelTime = timer.getTime() / num_meas;

	// format is: #groups #threads per group #copy time in mus #bandwidth in megabytes
	cout << ' ' << groups * threads_per_group << ' ' << groups << ' ' << threads_per_group << ' ' << kernelTime << ' ' << (2 * elems * sizeof(T) / kernelTime) << endl;
}
