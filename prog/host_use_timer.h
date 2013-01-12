/** @file
 * Time measurement utilities
 */
#ifndef _USETIMERH_
#define _USETIMERH_

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include "globaldefs.h"
#include "klepsydra/klepsydra.hpp"

extern char * benchmark_id;

/**
 * A wrapper for the monotonic Klepsydra timer
 * that can add up the time of multiple measurements.
 *
 * @todo Make this a child of klepsydra::Monotonic
 * @todo Merge this functionality back into Klepsydra.
 */
class usetimer {
public:
	/**
	 * Default and only constructor.
	 */
	usetimer() : time_measurement(0), num_meas(0) { };
	/**
	 * Reset the current time measurement.
	 */
	void reset();
	/**
	 * Store passed time and reset.
	 */
	void getTimeAndReset();
	/**
	 * Add passed time to any previously measured time and reset.
	 */
	void add();
	void add(uint64_t incr);
	/**
	 * Reset the aggregated measurement information.
	 */
	void zero();
	/**
	 * Retrieve the aggregated measured time in microseconds (10^6s).
	 */
	uint64_t getTime();
	/**
	 * Retrieve the number of measurements performed.
	 */
	int getNumMeas();
private:
	/**
	 * The aggregated measured time.
	 */
	uint64_t time_measurement;
	/**
	 * Currrently running measurement (not included in time_measurement)
	 */
	klepsydra::Monotonic timer;
	/**
	 * The number of measurements aggregated so far
	 */
	int num_meas;
};

/**
 * Measure the execution time of a OpenCL-Kernel based on an event associated to this kernel call known to have finished.
 */
uint64_t get_kernel_exec_time(cl_event event);

/**
 * Measure the time of a OpenCL-Kernel from queuing to calculation start based on an event associated to this kernel call known to have finished.
 */
uint64_t get_kernel_overhead_time(cl_event event);

/**
 * Measure the time an OpenCL-Kernel waited to be submitted to the device based on an event associated to this kernel call known to have finished.
 */
uint64_t get_kernel_submit_overhead_time(cl_event event);

/**
 * Divides a uint64_t-object by an integer.
 * Returns 0 if the integer is 0.
 */
uint64_t divide(uint64_t a, int b);

/**
 * Divides two uint64_t-objects and returns the result in percent.
 */
float percent(uint64_t a, uint64_t b);

/**
 * Callback function that resets the passed klepsydra timer on event completion.
 *
 * For use with clSetEventCallback
 */
void resetTimerOnComplete(cl_event event, cl_int event_command_exec_status, void * data);

#endif
