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

struct LatticeExtents2
{
	LatticeExtents2(SpatialLatticeExtent nsIn, TemporalLatticeExtent ntIn);
	LatticeExtents2(latticeSize nIn);
	LatticeExtents2(latticeSize nxIn, latticeSize nyIn, latticeSize nzIn, latticeSize ntIn);
	LatticeExtents2();
	const LatticeExtent xExtent, yExtent, zExtent, tExtent;
	latticeSize getNs() const;
	latticeSize getNt() const;
	latticeSize getSpatialLatticeVolume() const;
	latticeSize getLatticeVolume() const;
};

struct LatticeExtents
{
	LatticeExtents() : ns(4), nt(4) {};
	LatticeExtents(const int nsIn, const int ntIn): ns(nsIn), nt(ntIn) {};
	LatticeExtents(const unsigned int nsIn, const unsigned int ntIn): ns(nsIn), nt(ntIn) {};
	const unsigned int ns;
	const unsigned int nt;
	unsigned int getNs() const;
	unsigned int getNt() const;
	unsigned int getLatticeVolume() const;
	unsigned int getSpatialLatticeVolume() const;
};

struct BasicLatticeIndex
{
	BasicLatticeIndex(const latticeCoordinate x, const latticeCoordinate y, const latticeCoordinate z, const latticeCoordinate t, const LatticeExtents2 lE);
	BasicLatticeIndex up(const Direction dir) const;
	BasicLatticeIndex down(const Direction dir) const;
	operator latticeSize() const;
	const LatticeCoordinate x,y,z,t;
	const latticeSize spatialIndex, globalIndex;
};



