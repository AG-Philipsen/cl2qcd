/** @file
 * Implementation of the hardware::ProfilingData class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "profiling_data.hpp"
#include "system.hpp"

hardware::ProfilingData& hardware::ProfilingData::operator+=(const cl_event& event)
{
	cl_ulong time_start, time_end;
	cl_int err = clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &time_start, NULL);
	if(err) {
		throw hardware::OpenclException(err, "clGetEventProfilingInfo", __FILE__, __LINE__);
	}
	err = clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &time_end, NULL);
	if(err) {
		throw hardware::OpenclException(err, "clGetEventProfilingInfo", __FILE__, __LINE__);
	}

	total_time += (time_end - time_start) / 1e3;
	++num_values;

	return (*this);
}
