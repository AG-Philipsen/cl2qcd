/** @file
 * Declaration of the hardware::Device class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#ifndef _HARDWARE_DEVICE_HPP_
#define _HARDWARE_DEVICE_HPP_

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include "../meta/inputparameters.hpp"

namespace hardware {

	namespace buffers {
		// forward declaration for friend relation
		class Buffer;
	}

	/**
	 * Representation of an OpenCL device we are working on.
	 *
	 * Allows the querying and manipulation of the system and
	 * its hardware.
	 */
	class Device {

	friend hardware::buffers::Buffer;

	public:
		/**
		 * Initialize the device.
		 *
		 * \param context context for the device
		 * \param device_id id of the device to initialize
		 * \param inputparams the input parameters of the application
		 */
		Device(cl_context, cl_device_id, const meta::Inputparameters&);

		~Device();

		/**
		 * Checks whether the device supports double precision.
		 */
		bool is_double_supported();

		// non-copyable
		Device& operator=(const Device&) = delete;
		Device(const Device&) = delete;
		Device() = delete;

	private:
		/**
		 * The OpenCL context to be used by this device.
		 */
		cl_context context;

		/**
		 * The OpenCL device id of the device
		 */
		cl_device_id device_id;

		/**
		 * The input parameters of the application.
		 */
		const meta::Inputparameters& params;

		/**
		 * The command queue used to perform operations on this device
		 */
		cl_command_queue command_queue;
	};
}

#endif /* _HARDWARE_DEVICE_HPP_ */
