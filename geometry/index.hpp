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

#pragma once

#include "latticeExtents.hpp"
#include "../common_header_files/globaldefs.h"

typedef unsigned int uint;
typedef unsigned int latticeCoordinate;
typedef unsigned int latticeIndex;
typedef unsigned int latticeSize;

enum Direction {XDIR, YDIR, ZDIR, TDIR};

struct Index
{
	Index(const unsigned int, const unsigned int, const unsigned int, const unsigned int, const LatticeExtents);
	operator uint() const;
	const Index up(const Direction direction) const;
	const Index down(const Direction direction) const;
	latticeCoordinate x,y,z,t;
	LatticeExtents latticeExtents;
	latticeIndex globalIndex;
};

struct LinkIndex : public Index
{
	LinkIndex (const Index indexIn, const Direction dirIn);
	operator uint() const;
	const LinkIndex up(const Direction dirIn) const;
	const LinkIndex down(const Direction dirIn) const;
	Direction dir;
	latticeIndex globalIndex;
};
