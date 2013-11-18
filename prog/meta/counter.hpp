/** @file
 * Definition of a counter
 *
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
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

#ifndef _META_COUNTER_
#define _META_COUNTER_

namespace meta {

	/**
	 * Generic counter
	 */
	class Counter {

	public:
		/**
		 * Create the counter, initialized to 0
		 */
		Counter() noexcept;

		/**
		 * Increment the counter
		 */
		Counter& operator+=(const unsigned&) noexcept;

		/**
		 * Increment the counter
		 */
		Counter& operator++() noexcept;

		/**
		 * Evaluation operator
		 */
		operator unsigned() const noexcept;

		/**
		 * Reset the counter
		 */
		void reset() noexcept;

	private:
		/**
		 * The actual value of the counter.
		 */
		unsigned value;
	};
}

#endif /* _META_COUNTER_ */
