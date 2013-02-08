/** @file
 * Implementation of the hardware::System class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "system.hpp"

#include <list>
#include <sstream>
#include "../logger.hpp"
#include <stdexcept>

static std::vector<hardware::Device*> filter_cpus(const std::vector<hardware::Device*>& devices);

hardware::System::System(const meta::Inputparameters& params, bool enable_profiling)
	: params(params)
{
	using namespace hardware;

	// in debug scenarios make the compiler dump the compile results
	if( logger.beDebug() ) {
		setenv("GPU_DUMP_DEVICE_KERNEL", "3", 0); // can be overriden from outside
		setenv("AMD_OCL_BUILD_OPTIONS_APPEND", "-save-temps", 0); // can be overriden from outside
	}

	//LZ: for now, stick to one platform without any further checks...
	cl_int err = clGetPlatformIDs(1, &platform, 0);
	if(err) {
		throw OpenclException(err, "clGetPlatformIDs", __FILE__, __LINE__);
	}
	cl_context_properties context_props[3] = {
		CL_CONTEXT_PLATFORM,
		(cl_context_properties)platform,
		0
	};

	// restrict devices according to input parameters
	cl_device_type enabled_types = 0;
	if(params.get_use_gpu()) {
		enabled_types |= CL_DEVICE_TYPE_GPU;
	}
	if(params.get_use_cpu()) {
		enabled_types |= CL_DEVICE_TYPE_CPU;
	}

	// create devices
	context = clCreateContextFromType(context_props, enabled_types, 0, 0, &err);
	if(err) {
		throw OpenclException(err, "clCreateContextFromType", __FILE__, __LINE__);
	}
	cl_uint num_devices;
	err = clGetContextInfo(context, CL_CONTEXT_NUM_DEVICES, sizeof(cl_uint), &num_devices, 0);
	if(err) {
		throw OpenclException(err, "clGetContextInfo", __FILE__, __LINE__);
	}
	logger.info() << "Found " << num_devices << " OpenCL devices.";
	cl_device_id * device_ids = new cl_device_id[num_devices];
	err = clGetContextInfo(context, CL_CONTEXT_DEVICES, sizeof(cl_device_id) * num_devices, device_ids, 0);
	if(err) {
		throw OpenclException(err, "clGetContextInfo", __FILE__, __LINE__);
	}

	// check whether the user requested certain devices
	auto selection = params.get_selected_devices();
	if(selection.empty()) {
		// use all
		for(cl_uint i = 0; i < num_devices; ++i) {
			Device * dev = new Device(context, device_ids[i], params, enable_profiling);
#ifdef _USEDOUBLEPREC_
			if(!dev->is_double_supported()) {
				continue;
			}
#endif
			devices.push_back(dev);
		}
		// for now, if a gpu was found then throw out cpus
for(auto device: devices) {
			if(device->get_device_type() == CL_DEVICE_TYPE_GPU) {
				devices = filter_cpus(devices);
				break;
			}
		}
	} else {
for(int i: selection) {
			Device * dev = new Device(context, device_ids[i], params, enable_profiling);
#ifdef _USEDOUBLEPREC_
			if(!dev->is_double_supported()) {
				throw std::invalid_argument("Selected device does not support double precision.");
			}
#endif
			devices.push_back(dev);
		}
	}

	delete[] device_ids;
}

hardware::System::~System()
{
for(Device * device : devices) {
		delete device;
	}
	devices.clear();

	clReleaseContext(context);
}

const std::vector<hardware::Device*>& hardware::System::get_devices() const noexcept
{
	return devices;
}

const meta::Inputparameters& hardware::System::get_inputparameters() const noexcept
{
	return params;
}

hardware::OpenclException::OpenclException(int err)
{
	std::stringstream msg;
	msg << "OpenCL reported an error, error code: " << err;
	error_message = msg.str();
	logger.error() << error_message;
}

hardware::OpenclException::OpenclException(int err, std::string clname)
{
	std::stringstream msg;
	msg << "OpenCL reported an error in " << clname << ", error code: " << err;
	error_message = msg.str();
	logger.error() << error_message;
}

hardware::OpenclException::OpenclException(int err, std::string clname, std::string filename, int linenumber)
{
	std::stringstream msg;
	msg << "OpenCL failed. Error code " << err << " in " << clname << " at " << filename << ":" << linenumber;
	error_message = msg.str();
	logger.error() << error_message;
}

std::string hardware::OpenclException::what()
{
	return error_message;
}

hardware::System::operator const cl_context&() const noexcept
{
	return context;
}

void hardware::print_profiling(const System& system, const std::string& filename)
{
	auto devices = system.get_devices();
	for(size_t i = 0; i < devices.size(); ++i) {
		print_profiling(devices[i], filename, i);
	}
}

void hardware::print_profiling(const System * system, const std::string& filename)
{
	print_profiling(*system, filename);
}

static std::vector<hardware::Device*> filter_cpus(const std::vector<hardware::Device*>& devices)
{
	std::vector<hardware::Device*> filtered;
for(auto device: devices) {
		if(device->get_device_type() != CL_DEVICE_TYPE_CPU) {
			filtered.push_back(device);
		}
	}
	return filtered;
}
