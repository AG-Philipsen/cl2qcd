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
#define BOOST_TEST_MODULE geometry::Parallelization
#include <boost/test/unit_test.hpp>

#include "parallelization.hpp"

BOOST_AUTO_TEST_SUITE(TemporalParallelization)

	uint haloSize = 2;
	uint numberOfDevices = 4;
	uint ns = 4;
	uint nt = 8;
	LatticeExtents lE(ns, nt);
	LatticeGrid lG(numberOfDevices, lE);
	LocalLatticeExtents lLE(lE, lG);
	LatticeGridIndex lGI(0,0,0,0, lG);
	uint latticeGridPosition = 3;

	size_t sizePerElement = 15;

	size_t hyperVolumeLink = ns*ns*ns*4;
	size_t hyperVolumeNonLink = ns*ns*ns;
	size_t localExtentInSlowestDirection = nt/numberOfDevices;

	TemporalParallelizationHandlerLink tPHLink(lGI, lLE, sizePerElement, haloSize);
	TemporalParallelizationHandlerNonLink tPHNonLink(lGI, lLE, sizePerElement, haloSize);

	BOOST_AUTO_TEST_CASE(MAINPARTINDEX)
	{
		BOOST_REQUIRE_EQUAL( tPHLink.getMainPartIndex_destination(), 0);
		BOOST_REQUIRE_EQUAL( tPHNonLink.getMainPartIndex_destination(), 0);
	}

	BOOST_AUTO_TEST_CASE(FIRSTHALOPARTINDEX)
	{
		BOOST_REQUIRE_EQUAL( tPHLink.getFirstHaloIndex_destination(), hyperVolumeLink*localExtentInSlowestDirection);
		BOOST_REQUIRE_EQUAL( tPHNonLink.getFirstHaloIndex_destination(), hyperVolumeNonLink*localExtentInSlowestDirection);
	}

	BOOST_AUTO_TEST_CASE(SECONDHALOPARTINDEX)
	{
		BOOST_REQUIRE_EQUAL( tPHLink.getSecondHaloIndex_destination(), hyperVolumeLink*(localExtentInSlowestDirection + haloSize) );
		BOOST_REQUIRE_EQUAL( tPHNonLink.getSecondHaloIndex_destination(), hyperVolumeNonLink*(localExtentInSlowestDirection + haloSize) );
	}

	BOOST_AUTO_TEST_CASE(MAINPARTSIZE)
	{
		BOOST_REQUIRE_EQUAL( tPHLink.getMainPartSize(), hyperVolumeLink*localExtentInSlowestDirection );
		BOOST_REQUIRE_EQUAL( tPHNonLink.getMainPartSize(), hyperVolumeNonLink*localExtentInSlowestDirection );
	}

	BOOST_AUTO_TEST_CASE(HALOPARTSIZE)
	{
		BOOST_REQUIRE_EQUAL( tPHLink.getMainPartSize(), hyperVolumeLink*haloSize );
		BOOST_REQUIRE_EQUAL( tPHNonLink.getMainPartSize(), hyperVolumeNonLink*haloSize );
	}

	BOOST_AUTO_TEST_CASE(MAINPARTSIZEINBYTES)
	{
		BOOST_REQUIRE_EQUAL( tPHLink.getMainPartSizeInBytes(), hyperVolumeLink*localExtentInSlowestDirection*sizePerElement );
		BOOST_REQUIRE_EQUAL( tPHNonLink.getMainPartSizeInBytes(), hyperVolumeNonLink*localExtentInSlowestDirection*sizePerElement );
	}

	BOOST_AUTO_TEST_CASE(HALOPARTSIZEINBYTES)
	{
		BOOST_REQUIRE_EQUAL( tPHLink.getMainPartSizeInBytes(), hyperVolumeLink*haloSize*sizePerElement );
		BOOST_REQUIRE_EQUAL( tPHNonLink.getMainPartSizeInBytes(), hyperVolumeNonLink*haloSize*sizePerElement );
	}

	BOOST_AUTO_TEST_CASE(MAINPARTINDEX_SOURCE)
	{
		for(uint i = 0; i< numberOfDevices; i++)
		{
			LatticeGridIndex lGI(0,0,0,i, lG);
			TemporalParallelizationHandlerLink tPHLink(lGI, lLE, sizePerElement, haloSize);
			BOOST_REQUIRE_EQUAL( tPHLink.getMainPartIndex_source(), i* hyperVolumeLink*localExtentInSlowestDirection);
			TemporalParallelizationHandlerNonLink tPHNonLink(lGI, lLE, sizePerElement, haloSize);
			BOOST_REQUIRE_EQUAL( tPHNonLink.getMainPartIndex_source(), i* hyperVolumeNonLink*localExtentInSlowestDirection);
		}
	}

	BOOST_AUTO_TEST_CASE(FIRSTHALOPARTINDEX_SOURCE)
	{
		for(uint i = 0; i< numberOfDevices; i++)
		{
			LatticeGridIndex lGI(0,0,0,i, lG);
			TemporalParallelizationHandlerLink tPHLink(lGI, lLE, sizePerElement, haloSize);
			BOOST_REQUIRE_EQUAL( tPHLink.getFirstHaloPartIndex_source(), (i+1)%numberOfDevices* hyperVolumeLink*localExtentInSlowestDirection);
			TemporalParallelizationHandlerNonLink tPHNonLink(lGI, lLE, sizePerElement, haloSize);
			BOOST_REQUIRE_EQUAL( tPHNonLink.getFirstHaloPartIndex_source(), (i+1)%numberOfDevices* hyperVolumeNonLink*localExtentInSlowestDirection);
		}
	}

	BOOST_AUTO_TEST_CASE(SECONDHALOPARTINDEX_SOURCE)
	{
		for(uint i = 0; i< numberOfDevices; i++)
		{
			LatticeGridIndex lGI(0,0,0,i, lG);
			TemporalParallelizationHandlerLink tPHLink(lGI, lLE, sizePerElement, haloSize);
			BOOST_REQUIRE_EQUAL( tPHLink.getSecondHaloPartIndex_source(), hyperVolumeLink*((i*localExtentInSlowestDirection-haloSize+nt)%nt));
			TemporalParallelizationHandlerNonLink tPHNonLink(lGI, lLE, sizePerElement, haloSize);
			BOOST_REQUIRE_EQUAL( tPHNonLink.getSecondHaloPartIndex_source(), hyperVolumeNonLink*((i*localExtentInSlowestDirection-haloSize+nt)%nt));
		}
	}


BOOST_AUTO_TEST_SUITE_END()
