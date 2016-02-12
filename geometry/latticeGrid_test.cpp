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
#include <iostream>

#include "latticeGrid.hpp"

LatticeExtents lE(4, 8);
LatticeGrid lG(4);
LatticeGridIndex lGI(0,0,0,3,lG);
LocalLatticeExtents llE(lG, lE);
LocalLatticeMemoryExtents llME(lG, llE, 2);

BOOST_AUTO_TEST_SUITE(Grid)

	BOOST_AUTO_TEST_CASE(init)
	{
		BOOST_REQUIRE_EQUAL(lG.x, 1);
		BOOST_REQUIRE_EQUAL(lG.y, 1);
		BOOST_REQUIRE_EQUAL(lG.z, 1);
		BOOST_REQUIRE_EQUAL(lG.t, 4);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(GridNew)

	LatticeGridNew lG(4);

	BOOST_AUTO_TEST_CASE(init)
	{
		BOOST_REQUIRE_EQUAL(lG.xExtent, 1);
		BOOST_REQUIRE_EQUAL(lG.yExtent, 1);
		BOOST_REQUIRE_EQUAL(lG.zExtent, 1);
		BOOST_REQUIRE_EQUAL(lG.tExtent, 4);
	}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(GridIndex)

	BOOST_AUTO_TEST_CASE(exception)
	{
		BOOST_CHECK_THROW(LatticeGridIndex(1, 0, 0, 0, lG), std::logic_error);
		BOOST_CHECK_THROW(LatticeGridIndex(0, 1, 0, 0, lG), std::logic_error);
		BOOST_CHECK_THROW(LatticeGridIndex(0, 0, 1, 0, lG), std::logic_error);
		BOOST_CHECK_THROW(LatticeGridIndex(0, 0, 0, 4, lG), std::logic_error);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(LocalExtents)

	BOOST_AUTO_TEST_CASE(init)
	{
		BOOST_REQUIRE_EQUAL(llE.x, 4 / 1);
		BOOST_REQUIRE_EQUAL(llE.y, 4 / 1);
		BOOST_REQUIRE_EQUAL(llE.z, 4 / 1);
		BOOST_REQUIRE_EQUAL(llE.t, 8 / 4);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(LocalExtentsNew)

	uint numberOfDevices = 4;
	LatticeExtents2 tmp(SpatialLatticeExtent(4), TemporalLatticeExtent(8));
	LatticeGridNew lG(numberOfDevices);
	LocalLatticeExtentsNew llE(tmp, lG);

	BOOST_AUTO_TEST_CASE(init)
	{
		BOOST_REQUIRE_EQUAL(llE.xExtent, 4 / 1);
		BOOST_REQUIRE_EQUAL(llE.yExtent, 4 / 1);
		BOOST_REQUIRE_EQUAL(llE.zExtent, 4 / 1);
		BOOST_REQUIRE_EQUAL(llE.tExtent, 8 / numberOfDevices);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(LocalMemoryExtentsNew)

	uint numberOfDevices = 4;
	LatticeExtents2 tmp(SpatialLatticeExtent(4), TemporalLatticeExtent(8));
	LatticeGridNew lG(numberOfDevices);
	LocalLatticeExtentsNew llE(tmp, lG);
	uint haloSize = 2;
	LocalLatticeMemoryExtentsNew lME(lG, llE, haloSize);

	BOOST_AUTO_TEST_CASE(init)
	{
		BOOST_REQUIRE_EQUAL(lME.xExtent, llE.xExtent);
		BOOST_REQUIRE_EQUAL(lME.yExtent, llE.yExtent);
		BOOST_REQUIRE_EQUAL(lME.zExtent, llE.zExtent);
		BOOST_REQUIRE_EQUAL(lME.tExtent, llE.tExtent + 2*haloSize);
	}

	BOOST_AUTO_TEST_CASE(failure)
	{
		BOOST_REQUIRE_THROW(LocalLatticeMemoryExtentsNew(lG, llE, 20), std::invalid_argument);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(LocalMemoryExtents)

BOOST_AUTO_TEST_CASE(init)
{
	BOOST_REQUIRE_EQUAL(llME.x, 4);
	BOOST_REQUIRE_EQUAL(llME.y, 4);
	BOOST_REQUIRE_EQUAL(llME.z, 4);
	BOOST_REQUIRE_EQUAL(llME.t, 2 + 2 * 2);
}
BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(OverloadedStreamOperator)

BOOST_AUTO_TEST_CASE(LatticeGridStream)
{
	BOOST_CHECK_NO_THROW(std::cout << lG);
	BOOST_CHECK_NO_THROW(std::cout << lGI);
	BOOST_CHECK_NO_THROW(std::cout << llE);
	BOOST_CHECK_NO_THROW(std::cout << llME);
}
BOOST_AUTO_TEST_SUITE_END()
