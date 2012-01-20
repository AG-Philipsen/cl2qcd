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

void resetTimerOnComplete(cl_event event, cl_int event_command_exec_status, void * data)
{
	// complete is currently the only possible status, but just be sure
	if(event_command_exec_status == CL_COMPLETE) {
		klepsydra::Timer * timer = static_cast<klepsydra::Timer *>(data);
		timer->reset();
	}
}
