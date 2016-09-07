/*
 * Copyright 2016 Francesca Cuteri, Christopher Pinke
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

BOOST_AUTO_TEST_SUITE(LinkIndexBuild)

	BOOST_AUTO_TEST_CASE(linkIndex)
	{
		Direction dir = TDIR;
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( Index(1, 2, 3, 4, latticeExtents), dir)), NDIM * uint(Index(1, 2, 3, 4, latticeExtents)) + dir );
		dir = XDIR;
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( Index(1, 2, 3, 4, latticeExtents), dir)), NDIM * uint(Index(1, 2, 3, 4, latticeExtents)) + dir );
		dir = YDIR;
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( Index(1, 2, 3, 4, latticeExtents), dir)), NDIM * uint(Index(1, 2, 3, 4, latticeExtents)) + dir );
		dir = ZDIR;
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( Index(1, 2, 3, 4, latticeExtents), dir)), NDIM * uint(Index(1, 2, 3, 4, latticeExtents)) + dir );
	}

BOOST_AUTO_TEST_SUITE_END()


BOOST_AUTO_TEST_SUITE(LinkIndexNeighbours)

	Index index(3, 2, 1, 4, latticeExtents);

	BOOST_AUTO_TEST_CASE(upX)
	{
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, XDIR).up(XDIR)), XDIR + NDIM * (0 + 4*2 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, XDIR).up(YDIR)), XDIR + NDIM * (3 + 4*3 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, XDIR).up(ZDIR)), XDIR + NDIM * (3 + 4*2 + 4*4*2 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, XDIR).up(TDIR)), XDIR + NDIM * (3 + 4*2 + 4*4*1 + 4*4*4*0) );
	}

	BOOST_AUTO_TEST_CASE(upY)
	{
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, YDIR).up(XDIR)), YDIR + NDIM * (0 + 4*2 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, YDIR).up(YDIR)), YDIR + NDIM * (3 + 4*3 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, YDIR).up(ZDIR)), YDIR + NDIM * (3 + 4*2 + 4*4*2 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, YDIR).up(TDIR)), YDIR + NDIM * (3 + 4*2 + 4*4*1 + 4*4*4*0) );
	}

	BOOST_AUTO_TEST_CASE(upZ)
	{
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, ZDIR).up(XDIR)), ZDIR + NDIM * (0 + 4*2 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, ZDIR).up(YDIR)), ZDIR + NDIM * (3 + 4*3 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, ZDIR).up(ZDIR)), ZDIR + NDIM * (3 + 4*2 + 4*4*2 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, ZDIR).up(TDIR)), ZDIR + NDIM * (3 + 4*2 + 4*4*1 + 4*4*4*0) );
	}

	BOOST_AUTO_TEST_CASE(upT)
	{
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, TDIR).up(XDIR)), TDIR + NDIM * (0 + 4*2 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, TDIR).up(YDIR)), TDIR + NDIM * (3 + 4*3 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, TDIR).up(ZDIR)), TDIR + NDIM * (3 + 4*2 + 4*4*2 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, TDIR).up(TDIR)), TDIR + NDIM * (3 + 4*2 + 4*4*1 + 4*4*4*0) );
	}

	BOOST_AUTO_TEST_CASE(downX)
	{
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, XDIR).down(XDIR)), XDIR + NDIM * (2 + 4*2 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, XDIR).down(YDIR)), XDIR + NDIM * (3 + 4*1 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, XDIR).down(ZDIR)), XDIR + NDIM * (3 + 4*2 + 4*4*0 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, XDIR).down(TDIR)), XDIR + NDIM * (3 + 4*2 + 4*4*1 + 4*4*4*3) );
	}

	BOOST_AUTO_TEST_CASE(downY)
	{
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, YDIR).down(XDIR)), YDIR + NDIM * (2 + 4*2 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, YDIR).down(YDIR)), YDIR + NDIM * (3 + 4*1 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, YDIR).down(ZDIR)), YDIR + NDIM * (3 + 4*2 + 4*4*0 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, YDIR).down(TDIR)), YDIR + NDIM * (3 + 4*2 + 4*4*1 + 4*4*4*3) );
	}

	BOOST_AUTO_TEST_CASE(downZ)
	{
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, ZDIR).down(XDIR)), ZDIR + NDIM * (2 + 4*2 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, ZDIR).down(YDIR)), ZDIR + NDIM * (3 + 4*1 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, ZDIR).down(ZDIR)), ZDIR + NDIM * (3 + 4*2 + 4*4*0 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, ZDIR).down(TDIR)), ZDIR + NDIM * (3 + 4*2 + 4*4*1 + 4*4*4*3) );
	}

	BOOST_AUTO_TEST_CASE(downT)
	{
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, TDIR).down(XDIR)), TDIR + NDIM * (2 + 4*2 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, TDIR).down(YDIR)), TDIR + NDIM * (3 + 4*1 + 4*4*1 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, TDIR).down(ZDIR)), TDIR + NDIM * (3 + 4*2 + 4*4*0 + 4*4*4*4) );
		BOOST_REQUIRE_EQUAL(uint(LinkIndex( index, TDIR).down(TDIR)), TDIR + NDIM * (3 + 4*2 + 4*4*1 + 4*4*4*3) );
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(INDEX)

	BOOST_AUTO_TEST_CASE(SPATIAL_INDEX)
	{
		LatticeExtents lE(4,6,8,10);
		BOOST_REQUIRE_EQUAL( Index(1,2,3,4, lE).spatialIndex, 1 + 2*4 + 3*4*6 );
	}

	BOOST_AUTO_TEST_CASE(GLOBAL_INDEX)
	{
		LatticeExtents lE(4,6,8,10);
		BOOST_REQUIRE_EQUAL( Index(1,2,3,4, lE).globalIndex, 1 + 2*4 + 3*4*6 + 4*4*6*8);
	}

	uint extent = 4;
	LatticeExtents tmp(extent);

	BOOST_AUTO_TEST_CASE(NEIGHBOUR_UP_X)
	{
		for (uint i = 0; i< extent; i++)
		{
			BOOST_REQUIRE_EQUAL( Index(i,0,0,0, tmp).up(XDIR).x, (i+1)%extent );
		}
	}

	BOOST_AUTO_TEST_CASE(NEIGHBOUR_DOWN_X)
	{
		for (uint i = 0; i< extent; i++)
		{
			BOOST_REQUIRE_EQUAL( Index(i,0,0,0, tmp).down(XDIR).x, (i+extent-1)%extent );
		}
	}

	BOOST_AUTO_TEST_CASE(NEIGHBOUR_UP_Y)
	{
		for (uint i = 0; i< extent; i++)
		{
			BOOST_REQUIRE_EQUAL( Index(0,i,0,0, tmp).up(YDIR).y, (i+1)%extent );
		}
	}

	BOOST_AUTO_TEST_CASE(NEIGHBOUR_DOWN_Y)
	{
		for (uint i = 0; i< extent; i++)
		{
			BOOST_REQUIRE_EQUAL( Index(0,i,0,0, tmp).down(YDIR).y, (i+extent-1)%extent );
		}
	}

	BOOST_AUTO_TEST_CASE(NEIGHBOUR_UP_Z)
	{
		for (uint i = 0; i< extent; i++)
		{
			BOOST_REQUIRE_EQUAL( Index(0,0,i,0, tmp).up(ZDIR).z, (i+1)%extent );
		}
	}

	BOOST_AUTO_TEST_CASE(NEIGHBOUR_DOWN_Z)
	{
		for (uint i = 0; i< extent; i++)
		{
			BOOST_REQUIRE_EQUAL( Index(0,0,i,0, tmp).down(ZDIR).z, (i+extent-1)%extent );
		}
	}

	BOOST_AUTO_TEST_CASE(NEIGHBOUR_UP_T)
	{
		for (uint i = 0; i< extent; i++)
		{
			BOOST_REQUIRE_EQUAL( Index(0,0,0,i, tmp).up(TDIR).t, (i+1)%extent );
		}
	}

	BOOST_AUTO_TEST_CASE(NEIGHBOUR_DOWN_T)
	{
		for (uint i = 0; i< extent; i++)
		{
			BOOST_REQUIRE_EQUAL( Index(0,0,0,i, tmp).down(TDIR).t, (i+extent-1)%extent );
		}
	}

BOOST_AUTO_TEST_SUITE_END()


