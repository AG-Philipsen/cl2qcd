/** @file
 * Implementation of the hardware::Device class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "device.hpp"
#include "system.hpp"
#include "../logger.hpp"

static std::string retrieve_device_name(cl_device_id device_id);
static bool retrieve_device_availability(cl_device_id device_id);
static size_t retrieve_preferred_local_thread_num(cl_device_id device_id);
static size_t retrieve_preferred_global_thread_num(cl_device_id device_id);
static size_t retrieve_num_compute_units(cl_device_id device_id);
static cl_device_type retrieve_device_type(cl_device_id device_id);
static bool retrieve_supports_double(cl_device_id device_id);

hardware::Device::Device(cl_context context, cl_device_id device_id, const meta::Inputparameters& params, bool enable_profiling)
	: context(context), device_id(device_id), params(params),
	  preferred_local_thread_num(retrieve_preferred_local_thread_num(device_id)),
	  preferred_global_thread_num(retrieve_preferred_global_thread_num(device_id)),
	  num_compute_units(::retrieve_num_compute_units(device_id)),
	  device_type(::retrieve_device_type(device_id)),
	  supports_double(::retrieve_supports_double(device_id)),
	  prefers_blocked_loops(device_type == CL_DEVICE_TYPE_CPU),
	  name(retrieve_device_name(device_id)),
	  profiling_enabled(enable_profiling)
{
	logger.debug() << "Initializing " << retrieve_device_name(device_id);
	bool available = retrieve_device_availability(device_id);
	if(!available) {
		logger.error() << "Device is not available!";
	}

	cl_int err;
	logger.debug() << context << ' ' << device_id;
	command_queue = clCreateCommandQueue(context, device_id, profiling_enabled ? CL_QUEUE_PROFILING_ENABLE : 0, &err);
	if(err) {
		throw OpenclException(err, "clCreateCommandQueue", __FILE__, __LINE__);
	}
}

hardware::Device::~Device()
{
	clFinish(command_queue);
	clReleaseCommandQueue(command_queue);
}

bool hardware::Device::is_double_supported()
{
	return supports_double;
}

bool hardware::Device::get_prefers_blocked_loops() const noexcept
{
	return prefers_blocked_loops;
}


static std::string retrieve_device_name(cl_device_id device_id)
{
	using namespace hardware;
	size_t bytes;
	cl_int err = clGetDeviceInfo(device_id, CL_DEVICE_NAME, 0, 0, &bytes);
	if(err) {
		throw OpenclException(err, "clGetDeviceInfo(CL_DEVICE_NAME)", __FILE__, __LINE__);
	}
	char * name = new char[bytes + 1];
	err = clGetDeviceInfo(device_id, CL_DEVICE_NAME, bytes, name, 0);
	name[bytes] = 0;
	std::string val(name);
	delete[] name;
	if(err) {
		throw OpenclException(err, "clGetDeviceInfo(CL_DEVICE_NAME)", __FILE__, __LINE__);
	}
	return val;
}

static bool retrieve_device_availability(cl_device_id device_id)
{
	using namespace hardware;
	cl_bool available;
	cl_int err = clGetDeviceInfo(device_id, CL_DEVICE_AVAILABLE, sizeof(cl_bool), &available, 0);
	if(err) {
		throw OpenclException(err, "clGetDeviceInfo(CL_DEVICE_AVAILABLE)", __FILE__, __LINE__);
	}
	return available;
}

size_t hardware::Device::get_preferred_local_thread_num() const noexcept
{
	return preferred_local_thread_num;
}

size_t hardware::Device::get_preferred_global_thread_num() const noexcept
{
	return preferred_global_thread_num;
}

size_t hardware::Device::get_num_compute_units() const noexcept
{
	return num_compute_units;
}

cl_device_type hardware::Device::get_device_type() const noexcept
{
	return device_type;
}

static size_t retrieve_preferred_local_thread_num(cl_device_id device_id)
{
	if(retrieve_device_type(device_id) == CL_DEVICE_TYPE_GPU) {
		return 128;
	} else {
		return 1;
	}
}

static size_t retrieve_preferred_global_thread_num(cl_device_id device_id)
{
	size_t min_thread_num = retrieve_preferred_local_thread_num(device_id)
	                        * retrieve_num_compute_units(device_id);
	if(retrieve_device_type(device_id) == CL_DEVICE_TYPE_GPU) {
		return 4 * min_thread_num;
	} else {
		return min_thread_num;
	}
}

static size_t retrieve_num_compute_units(cl_device_id device_id)
{
	size_t num_compute_units;
	cl_int err = clGetDeviceInfo(device_id, CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(size_t), &num_compute_units, 0);
	if(err) {
		throw hardware::OpenclException(err, "clGetDeviceInfo(MAX_COMPUTE_UNITS)", __FILE__, __LINE__);
	}
	return num_compute_units;
}

static cl_device_type retrieve_device_type(cl_device_id device_id)
{
	cl_device_type device_type;
	cl_int err = clGetDeviceInfo(device_id, CL_DEVICE_TYPE, sizeof(cl_device_type), &device_type, 0);
	if(err) {
		throw hardware::OpenclException(err, "clGetDeviceInfo(TYPE)", __FILE__, __LINE__);
	}
	return device_type;
}

static bool retrieve_supports_double(cl_device_id device_id)
{
	using namespace hardware;
//  only on OpenCL 1.2
//	cl_device_fp_config double_support;
//	cl_int err = clGetDeviceInfo(device_id, CL_DEVICE_DOUBLE_FP_CONFIG, sizeof(double_support), &double_support, 0);
//	if(err) {
//		throw OpenclException();
//	}
//	return (double_support & (CL_FP_FMA | CL_FP_ROUND_TO_NEAREST | CL_FP_ROUND_TO_ZERO | CL_FP_ROUND_TO_INF | CL_FP_INF_NAN | CL_FP_DENORM));

	// backwards compatible query
	cl_int err;
	size_t value_size;
	err = clGetDeviceInfo(device_id, CL_DEVICE_EXTENSIONS, 0, 0, &value_size);
	if(err) {
		throw OpenclException(err, "clGetDeviceInfo(EXTENSIONS)", __FILE__, __LINE__);
	}
	char* extensions_val = new char[value_size + 1];
	err = clGetDeviceInfo(device_id, CL_DEVICE_EXTENSIONS, value_size, extensions_val, 0);
	extensions_val[value_size] = 0;
	std::string extensions(extensions_val);
	delete[] extensions_val;
	if(err) {
		throw OpenclException(err, "clGetDeviceInfo(EXTENSIONS)", __FILE__, __LINE__);
	}
	return (extensions.find("cl_khr_fp64") != std::string::npos);
}

hardware::Device::operator cl_command_queue() const noexcept
{
	return command_queue;
}

TmpClKernel hardware::Device::create_kernel(const char * const kernel_name, std::string build_opts) const
{
	if(params.is_ocl_compiler_opt_disabled()) {
		build_opts +=  " -cl-opt-disable";
	}
	return TmpClKernel(kernel_name, build_opts, context, &device_id, 1);
}

void hardware::Device::enqueue_kernel(cl_kernel kernel) const
{
	enqueue_kernel(kernel, get_preferred_global_thread_num());
}

void hardware::Device::enqueue_kernel(cl_kernel kernel, size_t global_threads) const
{
	enqueue_kernel(kernel, global_threads, get_preferred_local_thread_num());
}

void hardware::Device::enqueue_kernel(cl_kernel kernel, size_t global_threads, size_t local_threads) const
{
	cl_int clerr_enqueue = clEnqueueNDRangeKernel(command_queue, kernel, 1, 0, &global_threads, &local_threads, 0, 0, NULL);

	if(clerr_enqueue != CL_SUCCESS) {
		logger.fatal() << "clEnqueueNDRangeKernel failed, aborting...";
		logger.fatal() << "Some more information:";
		logger.fatal() << "global_work_size: " << global_threads;
		logger.fatal() << "local_work_size:  " << local_threads;

		size_t bytesInKernelName;
		if(clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, 0, NULL, &bytesInKernelName) == CL_SUCCESS) {
			char * kernelName = new char[bytesInKernelName]; // additional space for terminating 0 byte
			if(clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, bytesInKernelName, kernelName, NULL) == CL_SUCCESS) {
				logger.fatal() << "Failed kernel: " << kernelName;
			} else {
				logger.error() << "Could not retrieve kernel name";
			}
			delete [] kernelName;
		} else {
			logger.error() << "Could not retrieve length of kernel name";
		}

		throw hardware::OpenclException(clerr_enqueue, "clEnqueueNDRangeKernel", __FILE__, __LINE__);
	}
}

static int get_alignment_badness(size_t bytes)
{
	// this is pretty generic for all GPUs
	return (bytes % 256) ? 1 : 0;
}

static int get_cypress_stride_badness(size_t bytes, size_t lanes)
{
	const size_t CRITICAL_STRIDE = 16 * 1024; // at this stride performance is worst
	const size_t CRITICAL_STRIDE_RANGE = 768; // width of the critical stride

	int badness = 0;
	for(size_t hops = 1; hops < lanes; ++hops) {
		size_t dist_to_critical = (bytes * hops) % CRITICAL_STRIDE;

		if(dist_to_critical >= (CRITICAL_STRIDE - CRITICAL_STRIDE_RANGE) && dist_to_critical <= CRITICAL_STRIDE_RANGE) {
			++badness;
		}
	}
	return badness;
}

size_t hardware::Device::recommend_stride(size_t elems, size_t type_size, size_t lane_count) const
{
	size_t MAX_ADD_STRIDE = 8 * 1024; // never add more than 8 KiB per lane
	if(name == std::string("Cypress") || name == std::string("Cayman")) {
		// apply advanced stride rules
		for(size_t stride = elems; stride <= elems + MAX_ADD_STRIDE / type_size; ++stride) {
			if(get_cypress_stride_badness(stride * type_size, lane_count) == 0) {
				return stride;
			}
		}
		throw OptimizationError();
	} else {
		// simply align to 256 Bytes
		for(size_t stride = elems; stride <= elems + MAX_ADD_STRIDE / type_size; ++stride) {
			if(get_alignment_badness(stride * type_size) == 0) {
				return stride;
			}
		}
		throw OptimizationError();
	}
}

cl_device_id hardware::Device::get_id() const noexcept
{
	return device_id;
}

bool hardware::Device::is_profiling_enabled() const noexcept
{
	return profiling_enabled;
}
