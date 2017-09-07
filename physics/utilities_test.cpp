/** @file
 * Tests for functions working with sources.
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

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE physics
#include <boost/test/unit_test.hpp>

#include "utilities.hpp"

BOOST_AUTO_TEST_CASE(utilities)
{
	BOOST_REQUIRE_EQUAL(physics::buildCheckpointName( "conf.", "", 5, 1000),"conf.01000");
	BOOST_REQUIRE_EQUAL(physics::buildCheckpointName( "conf.", "postfix", 6, 1000),"conf.001000postfix");
	BOOST_REQUIRE_EQUAL(physics::buildCheckpointName( "conf.", "", 6, -1),"conf.save");
}
