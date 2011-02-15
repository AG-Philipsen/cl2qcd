/*
   This library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Library General Public
   License version 2 as published by the Free Software Foundation.

   This library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Library General Public License for more details.

   You should have received a copy of the GNU Library General Public License
   along with this library; see the file COPYING.LIB.  If not, write to
   the Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
   Boston, MA 02110-1301, USA.
*/

#include "host_timer.h"


Timer::Timer() {
	// make sure the timer is initialized with the current value
	reset();
}

void Timer::reset() {
	// replace the start value with the current value of the
	// monotonic clock
	start = getTimestamp();
}

uint64_t Timer::getTime() {
	uint64_t now = getTimestamp();
	
	return getDifference( start, now );
}

uint64_t Timer::getTimeAndReset() {
	uint64_t now = getTimestamp();
	
	uint64_t diff = getDifference( start, now );

	start = now;

	return diff;
}

uint64_t Timer::getDifference( uint64_t start, uint64_t end) const {
	uint64_t mus;
	mus = end - start;
	
	return mus;
}

uint64_t Timer::getTimestamp() const
{
	struct timespec now;		
	clock_gettime( CLOCK_MONOTONIC, &now );

	long nanos = now.tv_nsec;
	time_t seconds = now.tv_sec;
	
	uint64_t mus = static_cast<uint64_t>( seconds ) * 1000 * 1000 + nanos / 1000;
	
	return mus;
	
}
