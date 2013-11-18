/** @file
 * Implementation of the hardware::ProfilingData class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
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
