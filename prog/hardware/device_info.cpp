/** @file
 * Implementation of the hardware::DeviceInfo class
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
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

#include "device_info.hpp"
#include "system.hpp"

static std::string retrieve_device_name(cl_device_id device_id);
static size_t retrieve_preferred_local_thread_num(cl_device_id device_id);
static size_t retrieve_preferred_global_thread_num(cl_device_id device_id);
static size_t retrieve_num_compute_units(cl_device_id device_id);
static cl_device_type retrieve_device_type(cl_device_id device_id);


hardware::DeviceInfo::DeviceInfo(const cl_device_id device_id)
	: device_id(device_id),
	  preferred_local_thread_num(retrieve_preferred_local_thread_num(device_id)),
	  preferred_global_thread_num(retrieve_preferred_global_thread_num(device_id)),
	  num_compute_units(::retrieve_num_compute_units(device_id)),
	  device_type(::retrieve_device_type(device_id)),
	  prefers_blocked_loops(device_type == CL_DEVICE_TYPE_CPU),
	  prefers_soa(device_type == CL_DEVICE_TYPE_GPU),
	  name(retrieve_device_name(device_id))
{

}

hardware::DeviceInfo::DeviceInfo(const DeviceInfo& other)
	: device_id(other.device_id),
	  preferred_local_thread_num(other.preferred_local_thread_num),
	  preferred_global_thread_num(other.preferred_global_thread_num),
	  num_compute_units(other.num_compute_units),
	  device_type(other.device_type),
	  prefers_blocked_loops(other.prefers_blocked_loops),
	  prefers_soa(other.prefers_soa),
	  name(other.name)
{

}

bool hardware::DeviceInfo::is_double_supported() const noexcept
{
	return check_extension("cl_khr_fp64");
}

bool hardware::DeviceInfo::get_prefers_blocked_loops() const noexcept
{
	return prefers_blocked_loops;
}

bool hardware::DeviceInfo::get_prefers_soa() const noexcept
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

size_t hardware::DeviceInfo::get_preferred_local_thread_num() const noexcept
{
	return preferred_local_thread_num;
}

size_t hardware::DeviceInfo::get_preferred_global_thread_num() const noexcept
{
	return preferred_global_thread_num;
}

size_t hardware::DeviceInfo::get_num_compute_units() const noexcept
{
	return num_compute_units;
}

cl_device_type hardware::DeviceInfo::get_device_type() const noexcept
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

cl_device_id hardware::DeviceInfo::get_id() const noexcept
{
	return device_id;
}

std::string hardware::DeviceInfo::get_name() const noexcept
{
	return name;
}

bool hardware::DeviceInfo::check_extension(std::string const extension) const
{
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
	return (extensions.find(extension) != std::string::npos);

}
