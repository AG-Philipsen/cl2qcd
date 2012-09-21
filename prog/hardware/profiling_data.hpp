/** @file
 * Declaration of the hardware::Profiling class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#ifndef _HARDWARE_PROFILING_DATA_
#define _HARDWARE_PROFILING_DATA_

#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

namespace hardware {

	/**
	 * A utility class to aggregate profiling information
	 */
	class ProfilingData {

		public:
			/**
			 * Create an empty profiling data object
			 */
			ProfilingData() : ProfilingData(0, 0) {};

			/**
			 * Create an empty profiling data object
			 */
			ProfilingData(const ProfilingData& other) : ProfilingData(other.total_time, other.num_values) {};

			/**
			 * Query the total aggregated time
			 *
			 * \return time in ms
			 */
			inline uint64_t get_total_time() const noexcept {
				return total_time;
			}

			/**
			 * Query the number of aggregated values
			 *
			 * \return the number of aggregated values
			 */
			inline uint64_t get_num_values() const noexcept {
				return num_values;
			}

			/**
			 * Add a measurement.
			 *
			 * \param event An OpenCL profiling event
			 */
			ProfilingData& operator+=(const cl_event& event);

		private:
			/**
			 * Utility constructor used by the other constructors.
			 */
			ProfilingData(const uint64_t& total_time, const uint64_t num_values)
				: total_time(total_time), num_values(num_values) {};

			/**
			 * The total time aggregated in mus
			 */
			uint64_t total_time;

			/**
			 * The number of aggregated values
			 */
			uint64_t num_values;
	};
}

#endif /* _HARDWARE_PROFILING_DATA_ */
