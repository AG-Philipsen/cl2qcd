/**
 * @file
 *
 * Interface for monotonic timers.
 *
 * Copyright 2011 Matthias Bach <marix@marix.org>
 *
 * This file is part of Klepsydra.
 *
 * Klepsydra is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Klepsydra is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Klepsydra.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef _MONOTONIC_H_
#define _MONOTONIC_H_

#include "timer.hpp"

namespace klepsydra
{

	/**
	 * A timer that measures time using the monotonic
	 * timer at micro-second accuracy (if the system permits)
	 */
	class Monotonic : public Timer
	{
		public:
			/**
			 * Default and only constructor.
			 * Immediately starts the measurement.
			 */
			Monotonic();

			/**
			 * Resets the timer to start at 0 again.
			 */
			void reset();

			/**
			 * Recieves the time passed since timer start/reset.
			 * The value returned is in microseconds.
			 */
			uint64_t getTime();

			/**
			 * Recieves the time passed sind timer start/reset
			 * and resets the timer.
			 * The value returned is in microseconds.
			 */
			uint64_t getTimeAndReset();

		private:
			/**
			 * Time of the timer start/reset.
			 * This is the offset that needs to be substracted from
			 * later measurements to get the time difference.
			 */
			uint64_t start;

			/**
			 * Calculates the difference in microseconds between the two events
			 */
			uint64_t getDifference( uint64_t start, uint64_t end ) const;
		

			/**
			 * Retrieves the current Timestamp in an OS-specific manner
			 */
			uint64_t getTimestamp() const;
	};
}

#endif // _MONOTONIC_H_

