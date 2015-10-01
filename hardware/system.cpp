/** @file
 * Implementation of the hardware::System class
 *
 * Copyright (c) 2012-2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
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

/**
 * @Todo: Refactor/Remove dependence on meta::Inputparameters
 * This class is only needed for few things:
 * parameters.get_use_gpu(), parameters.get_use_cpu(), get_selected_devices(), get_device_count(), get_split_cpu()
 * and in creating the devices (which should not take meta::Inputparameters as well!
	new hardware::Device(context, info.get_id(), grid_pos, grid_size, params, enable_profiling))
 */

#include "system.hpp"

#include <list>
#include <sstream>
#include "../host_functionality/logger.hpp"
#include <stdexcept>
#include "device.hpp"
#include "transfer/transfer.hpp"

static std::list<hardware::DeviceInfo> filter_cpus(const std::list<hardware::DeviceInfo>& devices);
static std::vector<hardware::Device*> init_devices(const std::list<hardware::DeviceInfo>& infos, cl_context context, size_4 grid_size, const meta::Inputparameters& params, bool enable_profiling);
static size_4 calculate_grid_size(size_t num_devices);
static void setDebugEnvironmentVariables();

hardware::System::System(const meta::Inputparameters& params, const bool enable_profiling)
	: params(params), grid_size(0, 0, 0, 0), transfer_links()
{
	setDebugEnvironmentVariables();
	initOpenCLPlatforms();
	initOpenCLContext();
	initOpenCLDevices(enable_profiling);
}

void hardware::System::initOpenCLPlatforms()
{
	logger.debug() << "Init OpenCL platform(s)...";
	cl_uint numberOfAvailablePlatforms = 0;
	cl_int err = clGetPlatformIDs(1, &platform, &numberOfAvailablePlatforms);
	if(err)
	{
		throw OpenclException(err, "clGetPlatformIDs", __FILE__, __LINE__);
	}
	if (numberOfAvailablePlatforms > 1)
	{
		logger.warn() << "Found " << numberOfAvailablePlatforms << " platforms, take first one...";
	}
	else
	{
		logger.info() << "Found OpenCL platform";
	}
	logger.debug() << "...done";
}

cl_device_type restrictDeviceTypes( const meta::Inputparameters & parameters)
{
	// todo: does this cover all cases?
	cl_device_type enabled_types = 0;
	if(parameters.get_use_gpu()) {
		enabled_types |= CL_DEVICE_TYPE_GPU;
	}
	if(parameters.get_use_cpu()) {
		enabled_types |= CL_DEVICE_TYPE_CPU;
	}
	return enabled_types;
}

void hardware::System::initOpenCLContext()
{
	logger.debug() << "Init OpenCL context...";
	cl_int err = CL_SUCCESS;
	cl_context_properties context_props[3] = { CL_CONTEXT_PLATFORM, (cl_context_properties)platform, 0 };
	cl_device_type enabled_types = restrictDeviceTypes( params );

	context = clCreateContextFromType(context_props, enabled_types, 0, 0, &err);
	if(err)
	{
		throw OpenclException(err, "clCreateContextFromType", __FILE__, __LINE__);
	}
	logger.debug() << "...done";
}

void hardware::System::initOpenCLDevices(const bool enable_profiling)
{
	logger.debug() << "Init OpenCL devices...";
	cl_int err = CL_SUCCESS;

	cl_uint num_devices;
	err = clGetContextInfo(context, CL_CONTEXT_NUM_DEVICES, sizeof(cl_uint), &num_devices, 0);
	if(err)
	{
		throw OpenclException(err, "clGetContextInfo", __FILE__, __LINE__);
	}
	logger.info() << "Found " << num_devices << " OpenCL devices.";

	cl_device_id * device_ids = new cl_device_id[num_devices];
	err = clGetContextInfo(context, CL_CONTEXT_DEVICES, sizeof(cl_device_id) * num_devices, device_ids, 0);
	if(err)
	{
		throw OpenclException(err, "clGetContextInfo", __FILE__, __LINE__);
	}

	// check whether the user requested certain devices
	auto selection = params.get_selected_devices();
	std::list<DeviceInfo> device_infos;
	if(selection.empty()) {
		// use all (or up to max)
		size_t max_devices = params.get_device_count();
		for(cl_uint i = 0; i < num_devices && (!max_devices || device_infos.size() < max_devices); ++i) {
			DeviceInfo dev(device_ids[i]);
#ifdef _USEDOUBLEPREC_
			if(!dev.is_double_supported()) {
				logger.fatal() << "double not supported on device " << dev.get_name();
				continue;
			}
#endif
			device_infos.push_back(dev);
		}
		// for now, if a gpu was found then throw out cpus
for(auto device: device_infos) {
			if(device.get_device_type() == CL_DEVICE_TYPE_GPU) {
				device_infos = filter_cpus(device_infos);
				break;
			}
		}

		// if we are on a CPU, the number of devices is not restricted and we have OpenCL 1.2 split the CPU into NUMA domains (primarily makes testing easier)
#ifdef CL_VERSION_1_2
		if(params.get_split_cpu() && device_infos.size() == 1 && device_infos.front().get_device_type() == CL_DEVICE_TYPE_CPU && !max_devices) {
			cl_device_id original_device = device_infos.front().get_id();
			cl_device_partition_property partition_props[] = { CL_DEVICE_PARTITION_BY_AFFINITY_DOMAIN, CL_DEVICE_AFFINITY_DOMAIN_NEXT_PARTITIONABLE, 0};
			cl_uint num_sub_devs;
			err = clCreateSubDevices(original_device, partition_props, 0, nullptr, &num_sub_devs);
			if(err) {
				throw OpenclException(err, "clCreateSubDevices", __FILE__, __LINE__);
			}
			cl_device_id * sub_devices = new cl_device_id[num_sub_devs];
			err = clCreateSubDevices(original_device, partition_props, num_sub_devs, sub_devices, &num_sub_devs);
			if(err) {
				throw OpenclException(err, "clCreateSubDevices", __FILE__, __LINE__);
			}
			device_infos.clear();
			for(cl_uint i = 0; i < num_sub_devs; ++i) {
				device_infos.push_back(DeviceInfo(sub_devices[i]));
			}
			delete[] sub_devices;
		}
#endif
	} else {
		for(int i: selection) {
			if(i < 0 || i > (int) num_devices) {
				throw std::invalid_argument("Selected device does not exist");
			}
			DeviceInfo dev(device_ids[i]);
#ifdef _USEDOUBLEPREC_
			if(!dev.is_double_supported()) {
				throw std::invalid_argument("Selected device does not support double precision.");
			}
#endif
			device_infos.push_back(dev);
		}
	}

	grid_size = calculate_grid_size(device_infos.size());
	logger.info() << "Device grid layout: " << grid_size;

	devices = init_devices(device_infos, context, grid_size, params, enable_profiling);

	delete[] device_ids;

	logger.debug() << "...done";
}


hardware::System::~System()
{
	transfer_links.clear();

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
	errorCode = err;
	msg << "OpenCL reported an error, error code: " << err;
	error_message = msg.str();
	logger.error() << error_message;
}

hardware::OpenclException::OpenclException(int err, std::string clname)
{
	errorCode = err;
	std::stringstream msg;
	msg << "OpenCL reported an error in " << clname << ", error code: " << err;
	error_message = msg.str();
	logger.error() << error_message;
}

hardware::OpenclException::OpenclException(int err, std::string clname, std::string filename, int linenumber)
{
	errorCode = err;
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

static std::list<hardware::DeviceInfo> filter_cpus(const std::list<hardware::DeviceInfo>& devices)
{
	std::list<hardware::DeviceInfo> filtered;
	for(auto device: devices) {
		if(device.get_device_type() != CL_DEVICE_TYPE_CPU) {
			filtered.push_back(device);
		}
	}
	return filtered;
}

static std::vector<hardware::Device*> init_devices(const std::list<hardware::DeviceInfo>& infos, cl_context context, size_4 grid_size, const meta::Inputparameters& params, bool enable_profiling)
{
	std::vector<hardware::Device *> devices;
	devices.reserve(infos.size());

	unsigned tpos = 0;
	for(auto const info: infos) {
		size_4 const grid_pos(0, 0, 0, tpos++);
		if(grid_pos.t >= grid_size.t) {
			throw std::logic_error("Failed to place devices on the grid.");
		}
		devices.push_back(new hardware::Device(context, info.get_id(), grid_pos, grid_size, params, enable_profiling));
	}

	return devices;
}

static size_4 calculate_grid_size(size_t num_devices)
{
	// for now only parallelize in t-direction
	return size_4(1, 1, 1, num_devices);
}

size_4 hardware::System::get_grid_size()
{
	return grid_size;
}

hardware::Transfer * hardware::System::get_transfer(size_t from, size_t to, unsigned id) const
{
	auto const link_id = std::make_tuple(from, to, id);
	auto & link = transfer_links[link_id];
	if(!link.get()) {
		link = hardware::create_transfer(devices[from], devices[to], *this);
	}
	logger.trace() << "Serving Transfer: " << from << " -> " << to;
	return link.get();
}

cl_platform_id hardware::System::get_platform() const
{
	return platform;
}

static void makeCompilerDumbCompileResults()
{
	const bool overwriteExistingSettings = false;
	setenv("GPU_DUMP_DEVICE_KERNEL", "3", overwriteExistingSettings);
	setenv("AMD_OCL_BUILD_OPTIONS_APPEND", "-save-temps", overwriteExistingSettings);
}

static void setDebugEnvironmentVariables()
{
	if( logger.beDebug() ) {
		makeCompilerDumbCompileResults();
	}
}

