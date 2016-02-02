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

#include "index.hpp"

LatticeExtents latticeExtents(4, 5);

BOOST_AUTO_TEST_SUITE(IndexBuild)

	BOOST_AUTO_TEST_CASE(globalIndex)
	{
		BOOST_REQUIRE_EQUAL(uint(Index(1, 2, 3, 4, latticeExtents)), 1 + 4*2 + 4*4*3 + 4*4*4*4);
	}

	BOOST_AUTO_TEST_CASE(exception_1)
	{
		BOOST_CHECK_THROW(Index(5, 2, 3, 4, latticeExtents), std::invalid_argument);
		BOOST_CHECK_THROW(Index(1, 5, 3, 4, latticeExtents), std::invalid_argument);
		BOOST_CHECK_THROW(Index(1, 2, 5, 4, latticeExtents), std::invalid_argument);
		BOOST_CHECK_THROW(Index(1, 2, 3, 5, latticeExtents), std::invalid_argument);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(IndexNeighbours)

	BOOST_AUTO_TEST_CASE(upX)
	{
		BOOST_REQUIRE_EQUAL(uint(Index(1, 2, 3, 4, latticeExtents).up(XDIR)), 2 + 4*2 + 4*4*3 + 4*4*4*4);
	}

	BOOST_AUTO_TEST_CASE(upY)
	{
		BOOST_REQUIRE_EQUAL(uint(Index(1, 2, 3, 4, latticeExtents).up(YDIR)), 1 + 4*3 + 4*4*3 + 4*4*4*4);
	}

	BOOST_AUTO_TEST_CASE(upZ)
	{
		BOOST_REQUIRE_EQUAL(uint(Index(1, 2, 3, 4, latticeExtents).up(ZDIR)), 1 + 4*2 + 4*4*0 + 4*4*4*4);
	}

	BOOST_AUTO_TEST_CASE(upT)
	{
		BOOST_REQUIRE_EQUAL(uint(Index(1, 2, 3, 4, latticeExtents).up(TDIR)), 1 + 4*2 + 4*4*3 + 4*4*4*1);
	}

	BOOST_AUTO_TEST_CASE(downX)
	{
		BOOST_REQUIRE_EQUAL(uint(Index(1, 2, 3, 4, latticeExtents).down(XDIR)), 0 + 4*2 + 4*4*3 + 4*4*4*4);
	}

	BOOST_AUTO_TEST_CASE(downY)
	{
		BOOST_REQUIRE_EQUAL(uint(Index(1, 2, 3, 4, latticeExtents).down(YDIR)), 1 + 4*1 + 4*4*3 + 4*4*4*4);
	}

	BOOST_AUTO_TEST_CASE(downZ)
	{
		BOOST_REQUIRE_EQUAL(uint(Index(1, 2, 3, 4, latticeExtents).down(ZDIR)), 1 + 4*2 + 4*4*2 + 4*4*4*4);
	}

	BOOST_AUTO_TEST_CASE(downT)
	{
		BOOST_REQUIRE_EQUAL(uint(Index(1, 2, 3, 4, latticeExtents).down(TDIR)), 1 + 4*2 + 4*4*3 + 4*4*4*3);
	}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(LinkIndexBuild)

	BOOST_AUTO_TEST_CASE(linkIndex)
	{
		Direction dir = TDIR;
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( Index(1, 2, 3, 4, latticeExtents), dir)), uint(Index(1, 2, 3, 4, latticeExtents)) + latticeExtents.getLatticeVolume() * dir );
		dir = XDIR;
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( Index(1, 2, 3, 4, latticeExtents), dir)), uint(Index(1, 2, 3, 4, latticeExtents)) + latticeExtents.getLatticeVolume() * dir );
		dir = YDIR;
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( Index(1, 2, 3, 4, latticeExtents), dir)), uint(Index(1, 2, 3, 4, latticeExtents)) + latticeExtents.getLatticeVolume() * dir );
		dir = ZDIR;
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( Index(1, 2, 3, 4, latticeExtents), dir)), uint(Index(1, 2, 3, 4, latticeExtents)) + latticeExtents.getLatticeVolume() * dir );
	}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(LinkIndexNeighbours)

	Index index(3, 2, 1, 4, latticeExtents);

	BOOST_AUTO_TEST_CASE(upX)
	{
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, XDIR).up(XDIR)), XDIR*latticeExtents.getLatticeVolume() + (0 + 4*2 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, XDIR).up(YDIR)), XDIR*latticeExtents.getLatticeVolume() + (3 + 4*3 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, XDIR).up(ZDIR)), XDIR*latticeExtents.getLatticeVolume() + (3 + 4*2 + 4*4*2 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, XDIR).up(TDIR)), XDIR*latticeExtents.getLatticeVolume() + (3 + 4*2 + 4*4*1 + 4*4*4*1) );
	}

	BOOST_AUTO_TEST_CASE(upY)
	{
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, YDIR).up(XDIR)), YDIR*latticeExtents.getLatticeVolume() + (0 + 4*2 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, YDIR).up(YDIR)), YDIR*latticeExtents.getLatticeVolume() + (3 + 4*3 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, YDIR).up(ZDIR)), YDIR*latticeExtents.getLatticeVolume() + (3 + 4*2 + 4*4*2 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, YDIR).up(TDIR)), YDIR*latticeExtents.getLatticeVolume() + (3 + 4*2 + 4*4*1 + 4*4*4*1) );
	}

	BOOST_AUTO_TEST_CASE(upZ)
	{
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, ZDIR).up(XDIR)), ZDIR*latticeExtents.getLatticeVolume() + (0 + 4*2 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, ZDIR).up(YDIR)), ZDIR*latticeExtents.getLatticeVolume() + (3 + 4*3 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, ZDIR).up(ZDIR)), ZDIR*latticeExtents.getLatticeVolume() + (3 + 4*2 + 4*4*2 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, ZDIR).up(TDIR)), ZDIR*latticeExtents.getLatticeVolume() + (3 + 4*2 + 4*4*1 + 4*4*4*1) );
	}

	BOOST_AUTO_TEST_CASE(upT)
	{
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, TDIR).up(XDIR)), TDIR*latticeExtents.getLatticeVolume() + (0 + 4*2 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, TDIR).up(YDIR)), TDIR*latticeExtents.getLatticeVolume() + (3 + 4*3 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, TDIR).up(ZDIR)), TDIR*latticeExtents.getLatticeVolume() + (3 + 4*2 + 4*4*2 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, TDIR).up(TDIR)), TDIR*latticeExtents.getLatticeVolume() + (3 + 4*2 + 4*4*1 + 4*4*4*1) );
	}

	BOOST_AUTO_TEST_CASE(downX)
	{
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, XDIR).down(XDIR)), XDIR*latticeExtents.getLatticeVolume() + (2 + 4*2 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, XDIR).down(YDIR)), XDIR*latticeExtents.getLatticeVolume() + (3 + 4*1 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, XDIR).down(ZDIR)), XDIR*latticeExtents.getLatticeVolume() + (3 + 4*2 + 4*4*0 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, XDIR).down(TDIR)), XDIR*latticeExtents.getLatticeVolume() + (3 + 4*2 + 4*4*1 + 4*4*4*3) );
	}

	BOOST_AUTO_TEST_CASE(downY)
	{
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, YDIR).down(XDIR)), YDIR*latticeExtents.getLatticeVolume() + (2 + 4*2 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, YDIR).down(YDIR)), YDIR*latticeExtents.getLatticeVolume() + (3 + 4*1 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, YDIR).down(ZDIR)), YDIR*latticeExtents.getLatticeVolume() + (3 + 4*2 + 4*4*0 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, YDIR).down(TDIR)), YDIR*latticeExtents.getLatticeVolume() + (3 + 4*2 + 4*4*1 + 4*4*4*3) );
	}

	BOOST_AUTO_TEST_CASE(downZ)
	{
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, ZDIR).down(XDIR)), ZDIR*latticeExtents.getLatticeVolume() + (2 + 4*2 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, ZDIR).down(YDIR)), ZDIR*latticeExtents.getLatticeVolume() + (3 + 4*1 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, ZDIR).down(ZDIR)), ZDIR*latticeExtents.getLatticeVolume() + (3 + 4*2 + 4*4*0 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, ZDIR).down(TDIR)), ZDIR*latticeExtents.getLatticeVolume() + (3 + 4*2 + 4*4*1 + 4*4*4*3) );
	}

	BOOST_AUTO_TEST_CASE(downT)
	{
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, TDIR).down(XDIR)), TDIR*latticeExtents.getLatticeVolume() + (2 + 4*2 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, TDIR).down(YDIR)), TDIR*latticeExtents.getLatticeVolume() + (3 + 4*1 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, TDIR).down(ZDIR)), TDIR*latticeExtents.getLatticeVolume() + (3 + 4*2 + 4*4*0 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, TDIR).down(TDIR)), TDIR*latticeExtents.getLatticeVolume() + (3 + 4*2 + 4*4*1 + 4*4*4*3) );
	}


BOOST_AUTO_TEST_SUITE_END()
