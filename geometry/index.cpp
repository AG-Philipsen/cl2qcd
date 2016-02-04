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

Index::Index(const unsigned int xIn, const unsigned int yIn, const unsigned int zIn, const unsigned int tIn, const LatticeExtents lEIn ):
t(tIn), x(xIn), y(yIn), z(zIn), latticeExtents(lEIn)
{
	if( t >= lEIn.nt || x >= lEIn.ns || y >= lEIn.ns || z >= lEIn.ns )
	{
		throw std::invalid_argument( "Lattice indices must be smaller than lattice extents!" );
	}
	spaceIndex = x + y * lEIn.ns + z * lEIn.ns * lEIn.ns; // previously get_nspace
	globalIndex = x + y * lEIn.ns + z * lEIn.ns * lEIn.ns + t * lEIn.ns * lEIn.ns * lEIn.ns; // previously get_global_pos
}

Index::Index(const size_4 coordIn, const LatticeExtents lEIn):
t(coordIn.t), x(coordIn.x), y(coordIn.y), z(coordIn.z), latticeExtents(lEIn)
{
	if( t >= lEIn.nt || x >= lEIn.ns || y >= lEIn.ns || z >= lEIn.ns )
	{
		throw std::invalid_argument( "Lattice indices must be smaller than lattice extents!" );
	}
	spaceIndex = x + y * lEIn.ns + z * lEIn.ns * lEIn.ns; // previously get_nspace
	globalIndex = x + y * lEIn.ns + z * lEIn.ns * lEIn.ns + t * lEIn.ns * lEIn.ns * lEIn.ns; // previously get_global_pos
}

Index::operator uint() const
{
	return globalIndex;
}

const Index Index::up(const Direction dir) const
{
	switch(dir)
	{
	case XDIR:
		return Index( (x+1)%latticeExtents.ns, y, z, t, latticeExtents);
		break;

	case YDIR:
		return Index( x, (y+1)%latticeExtents.ns, z, t, latticeExtents);
		break;

	case ZDIR:
		return Index( x, y, (z+1)%latticeExtents.ns, t, latticeExtents);
		break;

	case TDIR:
		return Index( x, y, z, (t+1)%latticeExtents.ns, latticeExtents);
		break;
	default:
		throw std::invalid_argument( "Lattice direction must be between 0 and 3!" );
		break;
	}
}

const Index Index::down(const Direction dir) const
{
	switch(dir)
	{
	case XDIR:
		return Index( (x-1)%latticeExtents.ns, y, z, t, latticeExtents);
		break;

	case YDIR:
		return Index( x, (y-1)%latticeExtents.ns, z, t, latticeExtents);
		break;

	case ZDIR:
		return Index( x, y, (z-1)%latticeExtents.ns, t, latticeExtents);
		break;

	case TDIR:
		return Index( x, y, z, (t-1)%latticeExtents.ns, latticeExtents);
		break;
	default:
		throw std::invalid_argument( "Lattice direction must be between 0 and 3!" );
		break;
	}
}

LinkIndex::LinkIndex (const Index indexIn, const Direction dirIn):
	Index(indexIn), dir(dirIn)
{
	globalIndex = dir + NDIM * Index::globalIndex;
//	globalIndex = dir * latticeExtents.getLatticeVolume() + Index::globalIndex;
}

LinkIndex::operator uint() const
{
	return globalIndex;
}

const LinkIndex LinkIndex::up(const Direction dirIn) const
{
	return LinkIndex( Index::up(dirIn), dir );
}

const LinkIndex LinkIndex::down(const Direction dirIn) const
{
	return LinkIndex( Index::down(dirIn), dir );
}

uint LinkIndex::get_su3_idx_ildg_format(const uint n, const uint m)
{
	return 2 * n + 2 * m * NC + 2 * NC * NC * globalIndex;
}






