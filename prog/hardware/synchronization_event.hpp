/** @file
 * Declaration of the hardware::SynchronizationEvent class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#ifndef _HARDWARE_SYNCHRONIZATION_EVENT_HPP_
#define _HARDWARE_SYNCHRONIZATION_EVENT_HPP_

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif


namespace hardware {

/**
 * A wrapper class for a synchronization event.
 *
 * Allows to check whether a certain task has finished.
 */
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
	SynchronizationEvent& operator=(const SynchronizationEvent&) = delete;
	SynchronizationEvent(const SynchronizationEvent&);
	SynchronizationEvent() = delete; // default constructor does not make sense

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

	private:

	/**
	 * Reference to the wrapped OpenCL event.
	 */
	const cl_event event;
};

}

#endif /* _HARDWARE_SYNCHRONIZATION_EVENT_HPP_ */
