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

struct LatticeSpatialExtent
{
	LatticeSpatialExtent() : n(4) {};
	LatticeSpatialExtent(const int nsIn): n(nsIn) {};
	LatticeSpatialExtent(const unsigned int nsIn): n(nsIn) {};
	const unsigned int n;
};

struct LatticeTemporalExtent
{
	LatticeTemporalExtent() : n(4) {};
	LatticeTemporalExtent(const int nsIn): n(nsIn) {};
	LatticeTemporalExtent(const unsigned int nsIn): n(nsIn) {};
	const unsigned int n;
};

struct LatticeSymmetricExtent
{
	LatticeSymmetricExtent(): n(4) {};
	LatticeSymmetricExtent(const int nIn): n(nIn) {};
	LatticeSymmetricExtent(const unsigned int nIn): n(nIn) {};
	const unsigned int n;
};

struct LatticeExtents2
{
	LatticeExtents2(LatticeSpatialExtent nsIn, LatticeTemporalExtent ntIn) : ns(nsIn.n), nt(ntIn.n)  {};
	LatticeExtents2(LatticeSymmetricExtent nIn): ns(nIn.n), nt(nIn.n) {};
	LatticeSpatialExtent ns;
	LatticeTemporalExtent nt;
	unsigned int getNs() const;
	unsigned int getNt() const;
	unsigned int getLatticeVolume() const;
	unsigned int getSpatialLatticeVolume() const;
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
