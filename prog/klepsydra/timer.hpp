/**
 * @file
 *
 * Generic interface for Klepsydra timers.
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

#ifndef _TIMER_H_
#define _TIMER_H_

#include <stdint.h>

namespace klepsydra {

	/**
	 * Generic interface for all timers implemented by Klepsydra.
	 *
	 * TODO not all timers will be able to return time in microseconds,
	 * will have to be adjusted to be able to return alternate values
	 * like clock ticks.
	 */
	class Timer
	{
		public:
			/**
			 * Resets the timer to start at 0 again.
			 */
			virtual void reset() = 0;

			/**
			 * Recieves the time passed since timer start/reset.
			 * The value returned is in microseconds.
			 */
			virtual uint64_t getTime() = 0;

			/**
			 * Recieves the time passed sind timer start/reset
			 * and resets the timer.
			 * The value returned is in microseconds.
			 *
			 * The interface provides a basic implementations, but
			 * implementations are expected to provide optimized versions.
			 */
			virtual uint64_t getTimeAndReset();
	};

}

#endif // _TIMER_H_

