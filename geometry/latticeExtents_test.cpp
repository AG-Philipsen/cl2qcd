/**
 * Copyright 2016 Christopher Pinke
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
#define BOOST_TEST_MODULE geometry::latticeExtents
#include <boost/test/unit_test.hpp>

#include "latticeExtents.hpp"

BOOST_AUTO_TEST_SUITE(LATTICEEXTENT)

	BOOST_AUTO_TEST_CASE(BUILD1)
	{
		BOOST_REQUIRE_THROW( LatticeExtent(0), std::invalid_argument );
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(LATTICECOORDINATE)

	BOOST_AUTO_TEST_CASE(BUILD1)
	{
		BOOST_REQUIRE_THROW( LatticeCoordinate(4, 4), std::invalid_argument );
	}

	BOOST_AUTO_TEST_CASE(NEIGHBOUR_UP)
	{
		uint extent = 4;
		for (uint i = 0; i< extent; i++)
		{
			BOOST_REQUIRE_EQUAL( LatticeCoordinate(i, extent).up(), (i+1)%extent );
		}
	}

	BOOST_AUTO_TEST_CASE(NEIGHBOUR_DOWN)
	{
		uint extent = 4;
		for (uint i = 0; i< extent; i++)
		{
			BOOST_REQUIRE_EQUAL( LatticeCoordinate(i, extent).down(), (i-1+extent)%extent );
		}
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(LATTICEEXTENTS)

	LatticeExtents tmp(SpatialLatticeExtent(4), TemporalLatticeExtent(8));
	latticeSize ns = 4, nt=8;

	BOOST_AUTO_TEST_CASE(NS)
	{
		BOOST_REQUIRE_EQUAL( tmp.getNs(), 4 );
	}

	BOOST_AUTO_TEST_CASE(NT)
	{
		BOOST_REQUIRE_EQUAL( tmp.getNt(), 8 );
	}

	BOOST_AUTO_TEST_CASE(SPATIALVOLUME)
	{
		BOOST_REQUIRE_EQUAL( tmp.getSpatialLatticeVolume(), ns*ns*ns );
	}

	BOOST_AUTO_TEST_CASE(FROMFOURNUMBERS)
	{
		LatticeExtents tmp(2,3,4,8);
		BOOST_REQUIRE_EQUAL( tmp.getSpatialLatticeVolume(), 2*3*4 );
		BOOST_REQUIRE_EQUAL( tmp.getNs(), 2 );
	}

	BOOST_AUTO_TEST_CASE(VOLUME)
	{
		BOOST_REQUIRE_EQUAL( tmp.getLatticeVolume(), ns*ns*ns*nt );
	}

BOOST_AUTO_TEST_SUITE_END()

