/**
 * @file
 *
 * Implementation of the monotonic timer
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


#include "monotonic.hpp"

#ifdef __APPLE__
#include <mach/mach_time.h>

// method to get monotonic mac time, inspired by 
// http://www.wand.net.nz/~smr26/wordpress/2009/01/19/monotonic-time-in-mac-os-x/

#else
#include <time.h>
#endif


using namespace klepsydra;


Monotonic::Monotonic() {
	// make sure the timer is initialized with the current value
	reset();
}

void Monotonic::reset() {
	// replace the start value with the current value of the
	// monotonic clock
	start = getTimestamp();
}

uint64_t Monotonic::getTime() {
	uint64_t now = getTimestamp();

	return getDifference( start, now );
}

uint64_t Monotonic::getTimeAndReset() {
	uint64_t now = getTimestamp();

	uint64_t diff = getDifference( start, now );

	start = now;

	return diff;
}

uint64_t Monotonic::getDifference( uint64_t start, uint64_t end) const {
	uint64_t mus;
#ifdef __APPLE__
	uint64_t difference = end - start;
	static mach_timebase_info_data_t info = {0,0};

	if (info.denom == 0)
		mach_timebase_info(&info);

	uint64_t nanos = difference * (info.numer / info.denom);
	mus = nanos / 1000;
#else
	mus = end - start;
#endif

	return mus;
}

uint64_t Monotonic::getTimestamp() const
{
#ifdef __APPLE__
	return mach_absolute_time();
#else
	struct timespec now;
	clock_gettime( CLOCK_MONOTONIC, &now );

	long nanos = now.tv_nsec;
	time_t seconds = now.tv_sec;

	uint64_t mus = static_cast<uint64_t>( seconds ) * 1000 * 1000 + nanos / 1000;

	return mus;
#endif
}

