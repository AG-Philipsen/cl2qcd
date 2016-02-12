/*
 * Copyright 2015 Christopher Pinke
 * 			 2016 Francesca Cuteri
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

typedef unsigned int uint;
typedef unsigned int latticeSize;
typedef unsigned int latticeCoordinate;

enum Direction {TDIR = 0, XDIR, YDIR, ZDIR};

#include <ostream>
#include "../hardware/size_4.hpp"

struct LatticeExtent
{
	LatticeExtent(const latticeSize value);
	operator latticeSize() const;
	const latticeSize value;
};

struct LatticeCoordinate
{
	LatticeCoordinate(const latticeCoordinate valueIn, const LatticeExtent lE);
	LatticeCoordinate up() const;
	LatticeCoordinate down() const;
	operator latticeSize() const;

	const latticeCoordinate value;
	const LatticeExtent extent;
};

struct TemporalLatticeExtent;

struct SpatialLatticeExtent : public LatticeExtent
{
	SpatialLatticeExtent(const latticeSize value) : LatticeExtent(value){}
	SpatialLatticeExtent(const TemporalLatticeExtent&) = delete;
	SpatialLatticeExtent(TemporalLatticeExtent&) = delete;
};

struct TemporalLatticeExtent : public LatticeExtent
{
	TemporalLatticeExtent(const latticeSize value) : LatticeExtent(value){}
	TemporalLatticeExtent(const SpatialLatticeExtent&) = delete;
	TemporalLatticeExtent(SpatialLatticeExtent&) = delete;
};

struct LatticeExtents
{
	LatticeExtents(SpatialLatticeExtent nsIn, TemporalLatticeExtent ntIn);
	LatticeExtents(latticeSize nIn);
	LatticeExtents(latticeSize nxIn, latticeSize nyIn, latticeSize nzIn, latticeSize ntIn);
	LatticeExtents();
	const LatticeExtent xExtent, yExtent, zExtent, tExtent;
	latticeSize getNs() const;
	latticeSize getNt() const;
	latticeSize getSpatialLatticeVolume() const;
	latticeSize getLatticeVolume() const;
	operator size_4() const //todo: remove this!
		{
			return size_4(xExtent, yExtent, zExtent, tExtent);
		}
};

std::ostream& operator<<(std::ostream&, const LatticeExtents);
