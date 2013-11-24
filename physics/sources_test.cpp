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

#include "sources.hpp"
#include <sstream>

void test_sources(std::string type, int num_sources)
{
	using namespace physics::lattices;

	std::stringstream tmp;
	tmp << "--num_sources=";
	tmp << num_sources;
	std::string n_sources_string = tmp.str();
	std::string sourcetype_string = std::string("--sourcetype=") + type;
	const char * _params[] = {"foo", n_sources_string.c_str(), sourcetype_string.c_str()};
	meta::Inputparameters params(3, _params);
	hardware::System system(params);
	physics::PRNG prng(system);

	auto sources = create_sources(system, prng, params.get_num_sources());

	BOOST_REQUIRE_EQUAL(params.get_num_sources(), static_cast<const int>(sources.size()));

	release_spinorfields(sources);
}

BOOST_AUTO_TEST_CASE(sources)
{
	test_sources("point", 15);
	test_sources("volume", 2);
	test_sources("timeslice", 3);
	test_sources("zslice", 1);
}
