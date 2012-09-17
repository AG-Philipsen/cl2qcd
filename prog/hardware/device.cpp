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

hardware::Device::Device(cl_context context, cl_device_id device_id, const meta::Inputparameters& params)
	: context(context), device_id(device_id), params(params)
{
	logger.debug() << "Initializing " << retrieve_device_name(device_id);
	bool available = retrieve_device_availability(device_id);
	if(!available) {
		logger.error() << "Device is not available!";
	}

	cl_int err;
	logger.debug() << context << ' ' << device_id;
	command_queue = clCreateCommandQueue(context, device_id, 0, &err);
	if(err) {
		throw OpenclException(err, "clCreateCommandQueue", __FILE__, __LINE__);
	}
}

hardware::Device::~Device()
{
	clReleaseCommandQueue(command_queue);
}

/**
 * Checks whether the device supports double precision.
 */
bool hardware::Device::is_double_supported()
{
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
		throw OpenclException(err, "clGetDeviceInfo", __FILE__, __LINE__);
	}
	char* extensions_val = new char[value_size + 1];
	err = clGetDeviceInfo(device_id, CL_DEVICE_EXTENSIONS, value_size, extensions_val, 0);
	extensions_val[value_size] = 0;
	std::string extensions(extensions_val);
	delete[] extensions_val;
	if(err) {
		throw OpenclException(err, "clGetDeviceInfo", __FILE__, __LINE__);
	}
	return (extensions.find("cl_khr_fp64") != std::string::npos);
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

