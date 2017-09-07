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

#include "latticeGrid.hpp"
#include <stdexcept>
#include "../executables/exceptions.h"
#include "../host_functionality/logger.hpp"

LatticeGrid::LatticeGrid(const uint numberOfDevices, const LatticeExtents lE) :
	LatticeExtents(LatticeExtents(SpatialLatticeExtent(1),TemporalLatticeExtent(numberOfDevices)) )
{
	if(xExtent != 1 || yExtent != 1 || zExtent != 1)
	{
		throw std::invalid_argument("Only the time-direction can be parallelized");
	}
	if( lE.tExtent % numberOfDevices != 0)
	{
		throw std::invalid_argument("Cannot distribute lattice equally on devices!");
	}
}

LatticeGridIndex::LatticeGridIndex(const latticeCoordinate x, const latticeCoordinate y, const latticeCoordinate z, const latticeCoordinate t, const LatticeGrid lG) :
	Index(x,y,z,t,lG) {}

LocalLatticeExtents::LocalLatticeExtents(const LatticeExtents lE, const LatticeGrid lG) :
	LatticeExtents(lE.xExtent / lG.xExtent, lE.yExtent / lG.yExtent, lE.zExtent / lG.zExtent, lE.tExtent / lG.tExtent)
{
	if (xExtent % 2 || yExtent % 2 || zExtent %2 || tExtent % 2)
	{
		logger.warn() << "Local lattice size is odd. This is known to cause problems!";
	}
}

LocalLatticeMemoryExtents::LocalLatticeMemoryExtents(const LatticeGrid lG, const LocalLatticeExtents llE, unsigned int halo_size) :
	LatticeExtents(
			llE.xExtent + (lG.xExtent > 1 ? 2 * halo_size : 0), llE.yExtent + (lG.yExtent > 1 ? 2 * halo_size : 0),
			llE.zExtent + (lG.zExtent > 1 ? 2 * halo_size : 0), llE.tExtent + (lG.tExtent > 1 ? 2 * halo_size : 0) )
{
	if(llE.tExtent < halo_size) {
		throw std::invalid_argument("The lattice cannot be distributed onto the given grid.");
	}
}
