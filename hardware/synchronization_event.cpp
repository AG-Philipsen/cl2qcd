/** @file
 * Implementation of the hardware::SynchronizationEvent class
 *
 * Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
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

#include "synchronization_event.hpp"
#include "system.hpp"

hardware::SynchronizationEvent::SynchronizationEvent(const cl_event& event)
	: event(event)
{
	if(event) {
		cl_int err = clRetainEvent(event);
		if(err) {
			throw OpenclException(err, "clRetainEvent", __FILE__, __LINE__);
		}
	}
}

hardware::SynchronizationEvent::SynchronizationEvent(const hardware::SynchronizationEvent& other)
	: SynchronizationEvent(other.event) { }

hardware::SynchronizationEvent::SynchronizationEvent()
	: event(0) { }

hardware::SynchronizationEvent& hardware::SynchronizationEvent::operator=(const hardware::SynchronizationEvent& other)
{
	if(event) {
		cl_int err = clReleaseEvent(event);
		if(err) {
			throw OpenclException(err, "clRetainEvent", __FILE__, __LINE__);
		}
	}

	event = other.event;

	if(event) {
		cl_int err = clRetainEvent(event);
		if(err) {
			throw OpenclException(err, "clRetainEvent", __FILE__, __LINE__);
		}
	}
	return *this;
}

hardware::SynchronizationEvent::~SynchronizationEvent()
{
	if(event) {
		cl_int err = clReleaseEvent(event);
		if(err) {
			throw OpenclException(err, "clRetainEvent", __FILE__, __LINE__);
		}
	}
}

bool hardware::SynchronizationEvent::is_finished() const
{
	if(event) {
		cl_int event_status;
		cl_int err = clGetEventInfo(event, CL_EVENT_COMMAND_EXECUTION_STATUS, sizeof(cl_int), &event_status, 0);
		if(err) {
			throw OpenclException(err, "clGetEventInfo(CL_EVENT_COMMAND_EXECUTION_STATUS)", __FILE__, __LINE__);
		} else if(event_status < 0) {
			throw OpenclException(event_status, "event execution status", __FILE__, __LINE__);
		} else {
			return event_status == CL_COMPLETE;
		}
	} else {
		return true;
	}
}

void hardware::SynchronizationEvent::wait() const
{
	if(event) {
		cl_int err = clWaitForEvents(1, &event);
		if(err) {
			throw OpenclException(err, "clWaitForEvents", __FILE__, __LINE__);
		};
	}
}

void hardware::wait(const std::vector<SynchronizationEvent>& events)
{
	const size_t num_events = events.size();
	std::vector<cl_event> cl_events;
	cl_events.reserve(num_events);
	for(size_t i = 0; i < num_events; ++i) {
		if(events[i].is_valid()) {
			cl_events.push_back(events[i].raw());
		}
	}
	cl_int err = clWaitForEvents(num_events, cl_events.data());
	if(err) {
		throw OpenclException(err, "clWaitForEvents", __FILE__, __LINE__);
	};
}

const cl_event& hardware::SynchronizationEvent::raw() const
{
	return event;
}

bool hardware::SynchronizationEvent::is_valid() const
{
	return event != 0;
}

std::vector<cl_event> hardware::get_raw_events(const std::vector<SynchronizationEvent>& events)
{
	std::vector<cl_event> result;
	result.reserve(events.size());
	for(auto const & event: events) {
		if(event.is_valid()) {
			result.push_back(event.raw());
		}
	}
	return result;
}
