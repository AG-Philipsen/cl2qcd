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

FourUnsignedInt::FourUnsignedInt(const unsigned int x, const unsigned int y, const unsigned int z, const unsigned int t):
x(x), y(y), z(z), t(t)
{}

std::ostream& operator<<(std::ostream& out, const FourUnsignedInt & integers)
{
	return out << '(' << integers.x << ", " << integers.y << ", " << integers.z << ", " << integers.t << ')';
}

LatticeGrid::LatticeGrid( const unsigned int numberOfDevices):
FourUnsignedInt(1,1,1,numberOfDevices)
{
	if(x != 1 || y != 1 || z != 1) { //by now useless, but needed later on
		throw Print_Error_Message("Only the time-direction can be parallelized");
	}
}

LatticeGridIndex::LatticeGridIndex(const unsigned int x, const unsigned int y, const unsigned int z, const unsigned int t, const LatticeGrid lG):
FourUnsignedInt(x, y, z, t)
{
	if( x >= lG.x || y >= lG.y || z >= lG.z || t >= lG.t ) {
	throw std::logic_error("Failed to place devices on the grid.");
	}
}

LocalLatticeExtents::LocalLatticeExtents(const LatticeGrid lG, const LatticeExtents lE ):
		FourUnsignedInt(lE.getNs() / lG.x, lE.getNs() / lG.y, lE.getNs() / lG.z, lE.getNt() / lG.t)
{
	if (x % 2 || y % 2 || z %2 || t % 2)
	{
		logger.warn() << "Local lattice size is odd. This is known to cause problems!";
	}
	if (x * lG.x != lE.getNs() || y * lG.y != lE.getNs() || z * lG.z != lE.getNs() || t * lG.t != lE.getNt())
	{
		throw std::invalid_argument("The lattice cannot be distributed onto the given grid.");
	}

}

LocalLatticeMemoryExtents::LocalLatticeMemoryExtents( const LatticeGrid lG, const LocalLatticeExtents llE, unsigned int halo_size ):
		FourUnsignedInt(llE.x + (lG.x > 1 ? 2 * halo_size : 0), llE.y + (lG.y > 1 ? 2 * halo_size : 0), llE.z + (lG.z > 1 ? 2 * halo_size : 0), llE.t + (lG.t > 1 ? 2 * halo_size : 0))
{
	if(t < halo_size) {
		throw std::invalid_argument("The lattice cannot be distributed onto the given grid.");
	}
}
