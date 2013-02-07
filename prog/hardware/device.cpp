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
	  prefers_soa(device_type == CL_DEVICE_TYPE_GPU),
	  name(retrieve_device_name(device_id)),
	  profiling_enabled(enable_profiling),
	  profiling_data(),
	  gaugefield_code(nullptr),
	  prng_code(nullptr),
	  spinor_code(nullptr),
	  spinor_staggered_code(nullptr),
	  fermion_code(nullptr),
	  fermion_staggered_code(nullptr),
	  gaugemomentum_code(nullptr),
	  molecular_dynamics_code(nullptr),
	  correlator_code(nullptr),
	  heatbath_code(nullptr),
	  kappa_code(nullptr),
	  buffer_code(nullptr)
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
	if(buffer_code) {
		delete buffer_code;
	}
	if(kappa_code) {
		delete kappa_code;
	}
	if(heatbath_code) {
		delete heatbath_code;
	}
	if(correlator_code) {
		delete correlator_code;
	}
	if(gaugemomentum_code) {
		delete gaugemomentum_code;
	}
	if(fermion_code) {
		delete fermion_code;
	}
	if(spinor_code) {
		delete spinor_code;
	}
	if(prng_code) {
		delete prng_code;
	}
	if(gaugefield_code) {
		delete gaugefield_code;
	}

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

bool hardware::Device::get_prefers_soa() const noexcept
{
	return prefers_soa;
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
cl_command_queue hardware::Device::get_queue() const noexcept
{
	return command_queue;
}

TmpClKernel hardware::Device::create_kernel(const char * const kernel_name, std::string build_opts) const
{
	if(params.is_ocl_compiler_opt_disabled()) {
		build_opts +=  " -cl-opt-disable";
	}
	return TmpClKernel(kernel_name, build_opts, context, device_id);
}

void hardware::Device::enqueue_kernel(cl_kernel kernel)
{
	enqueue_kernel(kernel, get_preferred_global_thread_num());
}

void hardware::Device::enqueue_kernel(cl_kernel kernel, size_t global_threads)
{
	enqueue_kernel(kernel, global_threads, get_preferred_local_thread_num());
}

void hardware::Device::enqueue_kernel(cl_kernel kernel, size_t global_threads, size_t local_threads)
{
	// setup profiling if required
	cl_event profiling_event;
	// we only want to pass the event if we are actually profiling
	// otherwise the API will write back data into a no longer valid object
	cl_event * const profiling_event_p = profiling_enabled ? &profiling_event : 0;

	if(logger.beDebug() ) {
		logger.trace() << "calling clEnqueueNDRangeKernel...";
		logger.trace() << "global_work_size: " << global_threads;
		logger.trace() << "local_work_size:  " << local_threads;

		size_t bytesInKernelName;
		if(clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, 0, NULL, &bytesInKernelName) == CL_SUCCESS) {
			char * kernelName = new char[bytesInKernelName]; // additional space for terminating 0 byte
			if(clGetKernelInfo(kernel, CL_KERNEL_FUNCTION_NAME, bytesInKernelName, kernelName, NULL) == CL_SUCCESS) {
				logger.trace() << "Kernel: " << kernelName;
			} else {
				logger.error() << "Could not retrieve kernel name";
			}
			delete [] kernelName;
		} else {
			logger.error() << "Could not retrieve length of kernel name";
		}
	}

	// queue kernel
	cl_int clerr = clEnqueueNDRangeKernel(command_queue, kernel, 1, 0, &global_threads, &local_threads, 0, 0, profiling_event_p);

	// check for errors
	if(clerr != CL_SUCCESS) {
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

		throw hardware::OpenclException(clerr, "clEnqueueNDRangeKernel", __FILE__, __LINE__);
	}

	// evaluate profiling if required
	if(profiling_enabled) {
		// collect the data of this kernel invocation
		clerr = clWaitForEvents(1, &profiling_event);
		if(clerr) {
			throw hardware::OpenclException(clerr, "clWaitForEvents", __FILE__, __LINE__);
		}

		profiling_data[kernel] += profiling_event;
	}
}

void hardware::Device::enqueue_marker(cl_event * event) const
{
	cl_int err = clEnqueueMarker(command_queue, event);
	if(err) {
		throw hardware::OpenclException(err, "clEnqueueMarker()", __FILE__, __LINE__);
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

	int badness = get_alignment_badness(bytes);
	for(size_t hops = 1; hops < lanes; ++hops) {
		size_t dist_to_critical = (bytes * hops) % CRITICAL_STRIDE;

		if(dist_to_critical >= (CRITICAL_STRIDE - CRITICAL_STRIDE_RANGE) || dist_to_critical <= CRITICAL_STRIDE_RANGE) {
			++badness;
		}
	}
	return badness;
}

size_t hardware::Device::recommend_stride(size_t elems, size_t type_size, size_t lane_count) const
{
	size_t MAX_ADD_STRIDE = 8 * 1024; // never add more than 8 KiB per lane
	if(name == std::string("Cypress") || name == std::string("Cayman")) {
		logger.debug() << "Using cypress stride";
		// apply advanced stride rules
		for(size_t stride = elems; stride <= elems + MAX_ADD_STRIDE / type_size; ++stride) {
			if(get_cypress_stride_badness(stride * type_size, lane_count) == 0) {
				logger.debug() << "Return stride of " << stride << " elements, which is " << stride * type_size << " bytes.";
				return stride;
			}
		}
		throw OptimizationError();
	} else {
		// simply align to 256 Bytes
		logger.debug() << "Using default stride";
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

void hardware::Device::flush() const
{
	cl_int err = clFlush(command_queue);
	if(err) {
		throw hardware::OpenclException(err, "Failed when waiting for OpenCL device to finish.", __FILE__, __LINE__);
	}
}

void hardware::Device::synchronize() const
{
	cl_int err = clFinish(command_queue);
	if(err) {
		throw hardware::OpenclException(err, "Failed when waiting for OpenCL device to finish.", __FILE__, __LINE__);
	}
}

std::string hardware::Device::get_name() const noexcept
{
	return name;
}

hardware::ProfilingData hardware::Device::get_profiling_data(const cl_kernel& kernel) noexcept {
	return profiling_data[kernel];
}

hardware::code::Gaugefield * hardware::Device::get_gaugefield_code()
{
	if(!gaugefield_code) {
		gaugefield_code = new hardware::code::Gaugefield(params, this);
	}
	return gaugefield_code;
}

hardware::code::PRNG * hardware::Device::get_prng_code()
{
	if(!prng_code) {
		prng_code = new hardware::code::PRNG(params, this);
	}
	return prng_code;
}

hardware::code::Spinors * hardware::Device::get_spinor_code()
{
	if(!spinor_code) {
		spinor_code = new hardware::code::Spinors(params, this);
	}
	return spinor_code;
}

hardware::code::Spinors_staggered * hardware::Device::get_spinor_staggered_code()
{
	if(!spinor_staggered_code) {
		spinor_staggered_code = new hardware::code::Spinors_staggered(params, this);
	}
	return spinor_staggered_code;
}

hardware::code::Fermions * hardware::Device::get_fermion_code()
{
	if(!fermion_code) {
		fermion_code = new hardware::code::Fermions(params, this);
	}
	return fermion_code;
}

hardware::code::Fermions_staggered * hardware::Device::get_fermion_staggered_code()
{
	if(!fermion_staggered_code) {
		fermion_staggered_code = new hardware::code::Fermions_staggered(params, this);
	}
	return fermion_staggered_code;
}

hardware::code::Gaugemomentum * hardware::Device::get_gaugemomentum_code()
{
	if(!gaugemomentum_code) {
		gaugemomentum_code = new hardware::code::Gaugemomentum(params, this);
	}
	return gaugemomentum_code;
}

hardware::code::Molecular_Dynamics * hardware::Device::get_molecular_dynamics_code()
{
	if(!molecular_dynamics_code) {
		molecular_dynamics_code = new hardware::code::Molecular_Dynamics(params, this);
	}
	return molecular_dynamics_code;
}

hardware::code::Correlator * hardware::Device::get_correlator_code()
{
	if(!correlator_code) {
		correlator_code = new hardware::code::Correlator(params, this);
	}
	return correlator_code;
}

hardware::code::Heatbath * hardware::Device::get_heatbath_code()
{
	if(!heatbath_code) {
		heatbath_code = new hardware::code::Heatbath(params, this);
	}
	return heatbath_code;
}

hardware::code::Kappa * hardware::Device::get_kappa_code()
{
	if(!kappa_code) {
		kappa_code = new hardware::code::Kappa(params, this);
	}
	return kappa_code;
}

hardware::code::Buffer * hardware::Device::get_buffer_code()
{
	if(!buffer_code) {
		buffer_code = new hardware::code::Buffer(params, this);
	}
	return buffer_code;
}

void hardware::print_profiling(Device * device, const std::string& filename, int id)
{
	if(device->kappa_code) {
		device->kappa_code->print_profiling(filename, id);
	}
	if(device->heatbath_code) {
		device->heatbath_code->print_profiling(filename, id);
	}
	if(device->correlator_code) {
		device->correlator_code->print_profiling(filename, id);
	}
	if(device->fermion_code) {
		device->fermion_code->print_profiling(filename, id);
	}
	if(device->spinor_code) {
		device->spinor_code->print_profiling(filename, id);
	}
	if(device->prng_code) {
		device->prng_code->print_profiling(filename, id);
	}
	if(device->gaugefield_code) {
		device->gaugefield_code->print_profiling(filename, id);
	}
	if(device->buffer_code) {
		device->buffer_code->print_profiling(filename, id);
	}
}
