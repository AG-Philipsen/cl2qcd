/*
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
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

#include "host_use_timer.h"

#include "logger.hpp"

void usetimer::reset()
{
	timer.reset();
	return;
}

void usetimer::getTimeAndReset()
{
	time_measurement = timer.getTimeAndReset();
	//!!
	num_meas = 0;
	return;
}

void usetimer::add()
{
	time_measurement +=  timer.getTimeAndReset();
	num_meas ++;
	return;
}

void usetimer::add(uint64_t incr)
{
	time_measurement +=  incr;
	num_meas ++;
	return;
}

void usetimer::zero()
{
	time_measurement = 0;
	//!!
	num_meas = 0;
	return;
}

uint64_t usetimer::getTime()
{
	return time_measurement;
}

int usetimer::getNumMeas()
{
	return num_meas;
}

uint64_t divide(uint64_t a, int b)
{
	if (b == 0) return 0.;
	return uint64_t ( ( (float) a ) / ((float) b) );
}

float percent(uint64_t a, uint64_t b)
{
	if(b == 0) return 0;
	return (( (float) a) / ( (float )b )) * 100;
}

uint64_t get_kernel_exec_time(cl_event event)
{
	uint64_t tmp;
	cl_ulong time_start;
	cl_ulong time_end;
	clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &time_start, NULL);
	clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &time_end, NULL);
	tmp = (time_end - time_start ) * 0.001;
//  cout << tmp << endl;
	return tmp;
}

uint64_t get_kernel_overhead_time(cl_event event)
{
	uint64_t tmp;
	cl_ulong time_start;
	cl_ulong time_queue;
	clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &time_start, NULL);
	clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_QUEUED, sizeof(cl_ulong), &time_queue, NULL);
	tmp = (time_start - time_queue ) * 0.001;

	return tmp;
}

uint64_t get_kernel_submit_overhead_time(cl_event event)
{
	uint64_t tmp;
	cl_ulong time_submit;
	cl_ulong time_queue;
	clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_SUBMIT, sizeof(cl_ulong), &time_submit, NULL);
	clGetEventProfilingInfo(event, CL_PROFILING_COMMAND_QUEUED, sizeof(cl_ulong), &time_queue, NULL);
	tmp = (time_queue - time_submit ) * 0.001;

	return tmp;
}

void resetTimerOnComplete(cl_event, cl_int event_command_exec_status, void * data)
{
	// complete is currently the only possible status, but just be sure
	if(event_command_exec_status == CL_COMPLETE) {
		klepsydra::Timer * timer = static_cast<klepsydra::Timer *>(data);
		timer->reset();
	}
}
