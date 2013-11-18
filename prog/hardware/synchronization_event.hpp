/** @file
 * Declaration of the hardware::SynchronizationEvent class
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

#ifndef _HARDWARE_SYNCHRONIZATION_EVENT_HPP_
#define _HARDWARE_SYNCHRONIZATION_EVENT_HPP_

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif
#include <vector>


namespace hardware {

/**
 * A wrapper class for a synchronization event.
 *
 * Allows to check whether a certain task has finished.
 */
class SynchronizationEvent;

/**
 * Wait for multiple synchronization events.
 */
void wait(const std::vector<SynchronizationEvent>& events);

/**
 * Get the raw version of multiple events (filters out invalids)
 */
std::vector<cl_event> get_raw_events(const std::vector<SynchronizationEvent>& events);

class SynchronizationEvent {

public:

	/**
	 * Constructor
	 *
	 * @param event The OpenCL event to watch.
	 *              Don't forget to release it!
	 *              Class will claim it herself.
	 */
	SynchronizationEvent(const cl_event&);

	/**
	 * Release the event when class is destructed.
	 */
	~SynchronizationEvent();

	// make sure the event is aquired when copying
	SynchronizationEvent& operator=(const SynchronizationEvent&);
	SynchronizationEvent(const SynchronizationEvent&);
	SynchronizationEvent();

	/**
	 * Check whether the event has finished.
	 *
	 * Throws an exception if the event went into an aborted state.
	 */
	bool is_finished() const;

	/**
	 * Wait for the event to finish.
	 *
	 * Throws an exception if the event went into an aborted state.
	 */
	void wait() const;

	/**
	 * Provide a reference to the raw cl_event
	 */
	const cl_event& raw() const;

	/**
	 * Check whether this is a valid or a dummy event.
	 */
	bool is_valid() const;

private:

	/**
	 * Reference to the wrapped OpenCL event.
	 */
	cl_event event;
};

}

#endif /* _HARDWARE_SYNCHRONIZATION_EVENT_HPP_ */
