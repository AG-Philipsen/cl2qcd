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

#pragma once

#include "latticeExtents.hpp"
#include "../common_header_files/globaldefs.h"
//todo: remove!
#include "../hardware/size_4.hpp"

typedef uint latticeIndex;

struct Index
{
	Index(const unsigned int, const unsigned int, const unsigned int, const unsigned int, const LatticeExtents);
	Index(const size_4, const LatticeExtents);
	operator uint() const;
	const Index up(const Direction direction) const;
	const Index down(const Direction direction) const;
	latticeCoordinate t,x,y,z; //these should be const
	LatticeExtents latticeExtents;
	latticeIndex globalIndex, spaceIndex;
};
//todo: rename to Index
struct BasicLatticeIndex
{
	BasicLatticeIndex(const latticeCoordinate x, const latticeCoordinate y, const latticeCoordinate z, const latticeCoordinate t, const LatticeExtents2 lE);
	BasicLatticeIndex up(const Direction dir) const;
	BasicLatticeIndex down(const Direction dir) const;
	operator latticeSize() const;
	const LatticeCoordinate x,y,z,t;
	const latticeIndex spatialIndex, globalIndex;
};

struct LinkIndex : public Index
{
	LinkIndex (const Index indexIn, const Direction dirIn);
	operator latticeSize() const;
	const LinkIndex up(const Direction dirIn) const;
	const LinkIndex down(const Direction dirIn) const;
	uint get_su3_idx_ildg_format(const uint n, const uint m); // ADD TESTS!!
	const Direction direction;
	const latticeIndex globalIndex;
};

