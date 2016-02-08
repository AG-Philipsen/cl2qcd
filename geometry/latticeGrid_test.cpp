/*
 * Copyright 2016 Francesca Cuteri
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
#define BOOST_TEST_MODULE hardware::Device
#include <boost/test/unit_test.hpp>

#include "latticeGrid.hpp"

LatticeExtents lE(4, 8);
LatticeGrid lG(4);
LocalLatticeExtents llE(lG, lE);
LocalLatticeMemoryExtents llME(lG, llE, 2);

BOOST_AUTO_TEST_SUITE(Grid)

	BOOST_AUTO_TEST_CASE(init)
	{
		BOOST_REQUIRE_EQUAL(lG.nx, 1);
		BOOST_REQUIRE_EQUAL(lG.ny, 1);
		BOOST_REQUIRE_EQUAL(lG.nz, 1);
		BOOST_REQUIRE_EQUAL(lG.nt, 4);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(GridIndex)

	BOOST_AUTO_TEST_CASE(exception)
	{
		BOOST_CHECK_THROW(LatticeGridIndex(2, 1, 1, 4, lG), std::logic_error);
		BOOST_CHECK_THROW(LatticeGridIndex(1, 2, 1, 4, lG), std::logic_error);
		BOOST_CHECK_THROW(LatticeGridIndex(1, 1, 2, 4, lG), std::logic_error);
		BOOST_CHECK_THROW(LatticeGridIndex(1, 1, 1, 5, lG), std::logic_error);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(LocalExtents)

	BOOST_AUTO_TEST_CASE(init)
	{
		BOOST_REQUIRE_EQUAL(llE.nx, 4 / 1);
		BOOST_REQUIRE_EQUAL(llE.ny, 4 / 1);
		BOOST_REQUIRE_EQUAL(llE.nz, 4 / 1);
		BOOST_REQUIRE_EQUAL(llE.nt, 8 / 4);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(LocalMemoryExtents)

BOOST_AUTO_TEST_CASE(init)
{
	BOOST_REQUIRE_EQUAL(llME.nx, 4);
	BOOST_REQUIRE_EQUAL(llME.ny, 4);
	BOOST_REQUIRE_EQUAL(llME.nz, 4);
	BOOST_REQUIRE_EQUAL(llME.nt, 2 + 2 * 2);
}
BOOST_AUTO_TEST_SUITE_END()
