/** @file
 * Implementation of the hardware::System class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "system.hpp"

#include <list>
#include <sstream>
#include "../logger.hpp"

hardware::System::System(const meta::Inputparameters& params)
	: params(params)
{
	using namespace hardware;

	// in debug scenarios make the compiler dump the compile results
	if( logger.beDebug() ) {
		setenv("GPU_DUMP_DEVICE_KERNEL", "3", 0); // can be overriden from outside
		setenv("AMD_OCL_BUILD_OPTIONS_APPEND", "-save-temps", 0); // can be overriden from outside
	}

	//LZ: for now, stick to one platform without any further checks...
	cl_platform_id platform;
	cl_int err = clGetPlatformIDs(1, &platform, 0);
	if(err) {
		throw OpenclException(err, "clGetPlatformIDs", __FILE__, __LINE__);
	}
	cl_context_properties context_props[3] = {
		CL_CONTEXT_PLATFORM,
		(cl_context_properties)platform,
		0
	};

	// TODO allow restriction of device type via input params
	context = clCreateContextFromType(context_props, CL_DEVICE_TYPE_ALL, 0, 0, &err);
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

	for(cl_uint i = 0; i < num_devices; ++i) {
		Device * dev = new Device(context, device_ids[i], params);
#ifdef _USEDOUBLEPREC_
		if(!dev->is_double_supported()) {
			continue;
		}
#endif
		devices.push_back(dev);
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
	error_message = "OpenCL reported an error, error code: " + err;
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
