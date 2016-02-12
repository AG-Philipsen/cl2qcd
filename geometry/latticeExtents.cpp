/**
 * Copyright 2015 Christopher Pinke
 *           2016 Francesca Cuteri
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

#include "latticeExtents.hpp"
#include <stdexcept>

LatticeExtents::LatticeExtents(SpatialLatticeExtent nsIn, TemporalLatticeExtent ntIn) :
	xExtent(nsIn), yExtent(nsIn), zExtent(nsIn), tExtent(ntIn)	{}

LatticeExtents::LatticeExtents(latticeSize nIn):
	xExtent(nIn), yExtent(nIn), zExtent(nIn), tExtent(nIn)	{}

LatticeExtents::LatticeExtents(latticeSize nxIn, latticeSize nyIn, latticeSize nzIn, latticeSize ntIn):
	xExtent(nxIn), yExtent(nyIn), zExtent(nzIn), tExtent(ntIn)	{}

LatticeExtents::LatticeExtents() :
	xExtent(4), yExtent(4), zExtent(4), tExtent(4)	{}

latticeSize LatticeExtents::getNs() const
{
	return xExtent;
}

latticeSize LatticeExtents::getNt() const
{
	return tExtent;
}

latticeSize LatticeExtents::getSpatialLatticeVolume() const
{
	return xExtent * yExtent * zExtent;
}

latticeSize LatticeExtents::getLatticeVolume() const
{
	return getSpatialLatticeVolume() * tExtent;
}

LatticeExtent::LatticeExtent(const latticeSize value) : value(value)
{
	if(value == 0)
	{
		throw(std::invalid_argument("Must have non-zero latticeExtent!"));
	}
}

LatticeExtent::operator latticeSize() const
{
	return value;
}

LatticeCoordinate::LatticeCoordinate(const latticeCoordinate valueIn, const LatticeExtent lE):
	value(valueIn), extent(lE)
{
	if(value >= extent)
	{
		throw(std::invalid_argument("LatticeCoordinate must be smaller than LatticeExtent!"));
	}
}

LatticeCoordinate LatticeCoordinate::up() const
{
	return LatticeCoordinate( (value+1)%extent, extent);
}

LatticeCoordinate LatticeCoordinate::down() const
{
	return LatticeCoordinate( (value-1+extent)%extent, extent);
}

LatticeCoordinate::operator latticeSize() const
{
	return value;
}

std::ostream& operator<<(std::ostream& out, const LatticeExtents in)
{
	return out << '(' << in.xExtent << ", " << in.yExtent << ", " << in.zExtent << ", " << in.tExtent << ')';
}
