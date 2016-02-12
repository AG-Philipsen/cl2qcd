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

#include "index.hpp"
#include <stdexcept>

LinkIndex::LinkIndex (const Index indexIn, const Direction dirIn):
	Index(indexIn), direction(dirIn), globalIndex(direction + NDIM * Index::globalIndex) {}

LinkIndex::operator uint() const
{
	return globalIndex;
}

const LinkIndex LinkIndex::up(const Direction dirIn) const
{
	return LinkIndex( Index::up(dirIn), direction );
}

const LinkIndex LinkIndex::down(const Direction dirIn) const
{
	return LinkIndex( Index::down(dirIn), direction );
}

uint LinkIndex::get_su3_idx_ildg_format(const uint n, const uint m)
{
	return 2 * n + 2 * m * NC + 2 * NC * NC * globalIndex;
}

Index::Index(const latticeCoordinate x, const latticeCoordinate y, const latticeCoordinate z, const latticeCoordinate t, const LatticeExtents lE):
	x(x, lE.xExtent), y(y, lE.yExtent), z(z, lE.zExtent), t(t, lE.tExtent),
	spatialIndex(x + y*lE.xExtent + z*lE.xExtent*lE.yExtent),
	globalIndex(spatialIndex + t*lE.getSpatialLatticeVolume()) {}

Index::Index(const size_4 in, const LatticeExtents lE):
	x(in.x, lE.xExtent), y(in.y, lE.yExtent), z(in.z, lE.zExtent), t(in.t, lE.tExtent),
	spatialIndex(x + y*lE.xExtent + z*lE.xExtent*lE.yExtent),
	globalIndex(spatialIndex + t*lE.getSpatialLatticeVolume()) {}


Index Index::up(const Direction dir) const
{
	LatticeExtents tmp(x.extent, y.extent, z.extent, t.extent);
	switch(dir)
	{
	case XDIR:
		return Index( x.up(), y, z, t, tmp); break;

	case YDIR:
		return Index( x, y.up(), z, t, tmp); break;

	case ZDIR:
		return Index( x, y, z.up(), t, tmp); break;

	case TDIR:
		return Index( x, y, z, t.up(), tmp); break;

	default:
		throw std::invalid_argument( "Lattice direction must be between 0 and 3!" );
		break;
	}
}

Index Index::down(const Direction dir) const
{
	LatticeExtents tmp(x.extent, y.extent, z.extent, t.extent);
	switch(dir)
	{
	case XDIR:
		return Index( x.down(), y, z, t, tmp); break;

	case YDIR:
		return Index( x, y.down(), z, t, tmp); break;

	case ZDIR:
		return Index( x, y, z.down(), t, tmp); break;

	case TDIR:
		return Index( x, y, z, t.down(), tmp); break;

	default:
		throw std::invalid_argument( "Lattice direction must be between 0 and 3!" );
		break;
	}
}

Index::operator latticeSize() const
{
	return globalIndex;
}







