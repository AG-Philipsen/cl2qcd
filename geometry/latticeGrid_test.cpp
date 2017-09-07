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

BOOST_AUTO_TEST_SUITE(Grid)

	LatticeGrid lG(4, lE);

	BOOST_AUTO_TEST_CASE(init)
	{
		BOOST_REQUIRE_EQUAL(lG.xExtent, 1);
		BOOST_REQUIRE_EQUAL(lG.yExtent, 1);
		BOOST_REQUIRE_EQUAL(lG.zExtent, 1);
		BOOST_REQUIRE_EQUAL(lG.tExtent, 4);
	}

	BOOST_AUTO_TEST_CASE(FAILURE)
	{
		BOOST_REQUIRE_THROW( LatticeGrid(4, LatticeExtents(6)), std::invalid_argument );
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(GridIndex)

	LatticeGrid lG(4, lE);

	BOOST_AUTO_TEST_CASE(exception)
	{
		BOOST_CHECK_THROW(LatticeGridIndex(1, 0, 0, 0, lG), std::logic_error);
		BOOST_CHECK_THROW(LatticeGridIndex(0, 1, 0, 0, lG), std::logic_error);
		BOOST_CHECK_THROW(LatticeGridIndex(0, 0, 1, 0, lG), std::logic_error);
		BOOST_CHECK_THROW(LatticeGridIndex(0, 0, 0, 4, lG), std::logic_error);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(LocalExtents)

	uint numberOfDevices = 4;
	LatticeExtents tmp(SpatialLatticeExtent(4), TemporalLatticeExtent(8));
	LatticeGrid lG(numberOfDevices, tmp);
	LocalLatticeExtents llE(tmp, lG);

	BOOST_AUTO_TEST_CASE(init)
	{
		BOOST_REQUIRE_EQUAL(llE.xExtent, 4 / 1);
		BOOST_REQUIRE_EQUAL(llE.yExtent, 4 / 1);
		BOOST_REQUIRE_EQUAL(llE.zExtent, 4 / 1);
		BOOST_REQUIRE_EQUAL(llE.tExtent, 8 / numberOfDevices);
	}

BOOST_AUTO_TEST_SUITE_END()

BOOST_AUTO_TEST_SUITE(LocalMemoryExtents)

	uint numberOfDevices = 4;
	LatticeExtents tmp(SpatialLatticeExtent(4), TemporalLatticeExtent(8));
	LatticeGrid lG(numberOfDevices, tmp);
	LocalLatticeExtents llE(tmp, lG);
	uint haloSize = 2;
	LocalLatticeMemoryExtents lME(lG, llE, haloSize);

	BOOST_AUTO_TEST_CASE(init)
	{
		BOOST_REQUIRE_EQUAL(lME.xExtent, llE.xExtent);
		BOOST_REQUIRE_EQUAL(lME.yExtent, llE.yExtent);
		BOOST_REQUIRE_EQUAL(lME.zExtent, llE.zExtent);
		BOOST_REQUIRE_EQUAL(lME.tExtent, llE.tExtent + 2*haloSize);
	}

	BOOST_AUTO_TEST_CASE(failure)
	{
		BOOST_REQUIRE_THROW(LocalLatticeMemoryExtents(lG, llE, 20), std::invalid_argument);
	}

BOOST_AUTO_TEST_SUITE_END()


