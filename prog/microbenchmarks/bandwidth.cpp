/** @file
 * A simple microbenchmark to check device bandwidth for different
 * numbers of threada and streaming.
 */

#include <string>
#include <stdexcept>

#include <boost/program_options.hpp>

#include "../host_random.h"
#include "../opencl_module.h"
#include "../gaugefield_hybrid.h"
#include "../logger.hpp"
#include "../exceptions.h"

namespace po = boost::program_options;

extern std::string const version;
std::string const version = "0.1";

Random rnd(15);

/**
 * Selector type for the base type of the copy operations.
 */
enum copyType {
	type_invalid,
	type_float,
	type_su3,
	type_su3SOA,
	type_spinor,
	type_spinorSOA,
	type_spinorLocal
};

size_t getTypeSize(copyType type);

class Device : public Opencl_Module {

private:
	inputparameters params;
	cl_kernel floatKernel;
	cl_kernel su3Kernel;
	cl_kernel su3SOAKernel;
	cl_kernel spinorKernel;
	cl_kernel spinorSOAKernel;
	cl_kernel spinorLocalKernel;

	template<typename T> void runKernel(size_t groups, cl_ulong threads_per_group, cl_ulong elems, cl_kernel kernel, cl_mem in, cl_mem out);

public:
	Device(cl_command_queue queue, inputparameters* params, int maxcomp, string double_ext) : Opencl_Module() {
		Opencl_Module::init(queue, 0, params, maxcomp, double_ext); /* init in body for proper this-pointer */
	};
	virtual void fill_kernels();
	virtual void clear_kernels();
	~Device() {
		finalize();
	};

	void runKernel(copyType copy_type, size_t groups, cl_ulong threads_per_group, cl_ulong elems, cl_mem in, cl_mem out);
};

class Dummyfield : public Gaugefield_hybrid {

public:
	Dummyfield(cl_device_type device_type, size_t maxMemSize) : Gaugefield_hybrid(), maxMemSize(maxMemSize) {
		init(1, device_type, &params);
	};

	virtual void init_tasks();
	virtual void finalize_opencl();

	void runKernel(copyType copy_type, size_t groups, cl_ulong threads_per_group, cl_ulong elems);

private:
	void verify(hmc_complex, hmc_complex);
	void fill_buffers();
	void clear_buffers();
	inputparameters params;
	cl_mem in, out;
	size_t maxMemSize;
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
	type_map["spinor"] = type_spinor;
	type_map["spinorSOA"] = type_spinorSOA;
	type_map["spinorLocal"] = type_spinorLocal;

	const copyType copy_type = type_map[vm["type"].as<std::string>()];
	if(!copy_type) {
		logger.error() << "Please select one of the following types: float(default), su3, su3SOA, spinor, spinorSOA, spinorLocal";
		return 1;
	}
	logger.info() << "Using " << vm["type"].as<std::string>() << " as load/store datatype";

	if(vm.count("single") && vm.count("elements")) {
		logger.error() << "You can either use one element per thread or define the global number of elements";
	}

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
			Dummyfield dev(CL_DEVICE_TYPE_GPU, max_elements * getTypeSize(copy_type));
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
			Dummyfield dev(CL_DEVICE_TYPE_GPU, groups * max_threads * getTypeSize(copy_type));
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
			Dummyfield dev(CL_DEVICE_TYPE_GPU, elems * getTypeSize(copy_type));
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
			Dummyfield dev(CL_DEVICE_TYPE_GPU, max_groups * threads * getTypeSize(copy_type));
			for(size_t groups = 1; groups <= max_groups; ++groups) {
				dev.runKernel(copy_type, groups, threads, groups * threads);
			}
		} else {
			const cl_ulong elems = vm["elements"].as<cl_ulong>();
			logger.info() << "Keeping number of elements fixed at " << elems;
			Dummyfield dev(CL_DEVICE_TYPE_GPU, elems * getTypeSize(copy_type));
			for(size_t groups = 1; groups <= max_groups; ++groups) {
				dev.runKernel(copy_type, groups, threads, elems);
			}
		}
	}

	logger.info() << "Done";

	return 0;
}

void Dummyfield::fill_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	cl_int err;

	cl_context context = opencl_modules[0]->get_context();

	logger.info() << "Allocating buffers of " << maxMemSize << " bytes.";

	in = clCreateBuffer(context, CL_MEM_READ_ONLY, maxMemSize, 0, &err );
	if(err) {
		logger.fatal() << "Unable to allocate memory on device";
		throw Opencl_Error(err);
	}

	out = clCreateBuffer(context, CL_MEM_WRITE_ONLY, maxMemSize, 0, &err );
	if(err) {
		logger.fatal() << "Unable to allocate memory on device";
		throw Opencl_Error(err);
	}
}

void Device::fill_kernels()
{
	Opencl_Module::fill_kernels();

	floatKernel = createKernel("copyFloat") << basic_opencl_code << "types_fermions.h" << "microbenchmarks/bandwidth.cl";
	su3Kernel = createKernel("copySU3") << basic_opencl_code << "types_fermions.h" << "microbenchmarks/bandwidth.cl";
	su3SOAKernel = createKernel("copySU3SOA") << basic_opencl_code << "types_fermions.h" << "microbenchmarks/bandwidth.cl";
	spinorKernel = createKernel("copySpinor") << basic_opencl_code << "types_fermions.h" << "microbenchmarks/bandwidth.cl";
	spinorSOAKernel = createKernel("copySpinorSOA") << basic_opencl_code << "types_fermions.h" << "microbenchmarks/bandwidth.cl";
	spinorLocalKernel = createKernel("copySpinorLocal") << basic_opencl_code << "types_fermions.h" << "microbenchmarks/bandwidth.cl";
}

void Dummyfield::clear_buffers()
{
	// don't invoke parent function as we don't require the original buffers

	clReleaseMemObject(in);
	clReleaseMemObject(out);
}

void Device::clear_kernels()
{
	// don't invoke parent function as we don't require the original kernels

	clReleaseKernel(floatKernel);
	clReleaseKernel(su3Kernel);
	clReleaseKernel(su3SOAKernel);
	clReleaseKernel(spinorKernel);
	clReleaseKernel(spinorSOAKernel);
	clReleaseKernel(spinorLocalKernel);
}


void Device::runKernel(copyType copy_type, size_t groups, cl_ulong threads_per_group, cl_ulong elems, cl_mem in, cl_mem out)
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
		case type_spinor:
			runKernel<spinor>(groups, threads_per_group, elems, spinorKernel, in, out);
			return;
		case type_spinorSOA:
			runKernel<spinor>(groups, threads_per_group, elems, spinorSOAKernel, in, out);
			return;
		case type_spinorLocal:
			runKernel<spinor>(groups, threads_per_group, elems, spinorLocalKernel, in, out);
			return;
		default:
			throw invalid_argument("runKernel has not been implemented for this type");
	}
}

template<typename T> void Device::runKernel(size_t groups, cl_ulong threads_per_group, cl_ulong elems, cl_kernel kernel, cl_mem in, cl_mem out)
{
	cl_int err = CL_SUCCESS;

	// TODO adjust threads_per_group for kernel invocation to proper size but keep requested value for kernel arg
	size_t local_threads = threads_per_group;
	size_t total_threads = groups * local_threads;

	err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &out);
	if(err) {
		logger.fatal() << "Failed to set kernel argument: " << err;
		throw Opencl_Error(err);
	}
	err = clSetKernelArg(kernel, 1, sizeof(cl_mem), &in);
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

	cl_command_queue queue = get_queue();

	err = clEnqueueNDRangeKernel(get_queue(), kernel, 1, 0, &total_threads, &local_threads, 0, 0, NULL);
	err = clFinish(get_queue());
	if(err) {
		logger.fatal() << "Failed to execute kernel: " << err;
		throw Opencl_Error(err);
	}
	klepsydra::Monotonic timer;
	for(size_t i = 0; i < num_meas; ++i)
		clEnqueueNDRangeKernel(get_queue(), kernel, 1, 0, &total_threads, &local_threads, 0, 0, NULL);
	err = clFinish(queue);
	int64_t kernelTime = timer.getTime() / num_meas;
	if(err) {
		logger.fatal() << "Failed to execute kernel: " << err;
		throw Opencl_Error(err);
	}

	// format is: #groups #threads per group #elements #copied memory in bytes #copy time in mus #bandwidth in megabytes
	// FIXME sizeof can give broken results in case of aligned types (gross size not equal to net content size)
	cout << groups * threads_per_group << ' ' << groups << ' ' << threads_per_group << ' ' << elems << ' ' << elems * sizeof(T) << ' ' << kernelTime << ' ' << (2 * elems * sizeof(T) / kernelTime) << endl;
}

void Dummyfield::init_tasks()
{
	opencl_modules = new Opencl_Module* [get_num_tasks()];
	opencl_modules[0] = new Device(queue[0], get_parameters(), get_max_compute_units(0), get_double_ext(0));

	fill_buffers();
}

void Dummyfield::finalize_opencl()
{
	clear_buffers();
	Gaugefield_hybrid::finalize_opencl();
}

void Dummyfield::runKernel(copyType copy_type, size_t groups, cl_ulong threads_per_group, cl_ulong elems)
{
	static_cast<Device*>(opencl_modules[0])->runKernel(copy_type, groups, threads_per_group, elems, in, out);
}

size_t getTypeSize(copyType type)
{
	switch(type) {
		case type_float:
			return sizeof(hmc_float);
		case type_su3:
		case type_su3SOA:
			return sizeof(Matrixsu3);
		case type_spinor:
		case type_spinorSOA:
		case type_spinorLocal:
			return sizeof(spinor);
		default:
			throw invalid_argument("getTypeSize has not been implemented for this type.");
	}
}
