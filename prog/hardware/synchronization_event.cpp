/** @file
 * Implementation of the hardware::SynchronizationEvent class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "synchronization_event.hpp"
#include "system.hpp"

hardware::SynchronizationEvent::SynchronizationEvent(const cl_event& event)
	: event(event)
{
	cl_int err = clRetainEvent(event);
	if(err) {
		throw OpenclException(err, "clRetainEvent", __FILE__, __LINE__);
	}
}

hardware::SynchronizationEvent::SynchronizationEvent(const hardware::SynchronizationEvent& other)
	: SynchronizationEvent(other.event) { }

hardware::SynchronizationEvent::~SynchronizationEvent()
{
	cl_int err = clReleaseEvent(event);
	if(err) {
		throw OpenclException(err, "clRetainEvent", __FILE__, __LINE__);
	}
}

bool hardware::SynchronizationEvent::is_finished() const
{
	cl_int event_status;
	cl_int err = clGetEventInfo(event, CL_EVENT_COMMAND_EXECUTION_STATUS, sizeof(cl_int), &event_status, 0);
	if(err) {
		throw OpenclException(err, "clGetEventInfo(CL_EVENT_COMMAND_EXECUTION_STATUS)", __FILE__, __LINE__);
	} else if(event_status < 0) {
		throw OpenclException(event_status, "event execution status", __FILE__, __LINE__);
	} else {
		return event_status == CL_COMPLETE;
	}
}

void hardware::SynchronizationEvent::wait() const
{
	cl_int err = clWaitForEvents(1, &event);
	if(err) {
		throw OpenclException(err, "clWaitForEvents", __FILE__, __LINE__);
	};
}
