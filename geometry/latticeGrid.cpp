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

LatticeGrid::LatticeGrid( const unsigned int numberOfDevices):
nx(1), ny(1), nz(1), nt(numberOfDevices)
{
	if(nx != 1 || ny != 1 || nz != 1) { //by now useless, but needed later on
		throw Print_Error_Message("Only the time-direction can be parallelized");
	}
}

LatticeGridIndex::LatticeGridIndex(const unsigned int x, const unsigned int y, const unsigned int z, const unsigned int t, const LatticeGrid lG):
x(x), y(y), z(z), t(t)
{
	if( x >= lG.nx || y >= lG.ny || z >= lG.nz || t >= lG.nt ) {
	throw std::logic_error("Failed to place devices on the grid.");
	}
}

LocalLatticeExtents::LocalLatticeExtents( LatticeGrid lG, LatticeExtents lE ):
		nx(lE.getNs() / lG.nx),
		ny(lE.getNs() / lG.ny),
		nz(lE.getNs() / lG.nz),
		nt(lE.getNt() / lG.nt)
{
	if (nx % 2 || ny % 2 || nz %2 || nt % 2)
	{
		logger.warn() << "Local lattice size is odd. This is known to cause problems!";
	}
	if (nx * lG.nx != lE.getNs() || ny * lG.ny != lE.getNs() || nz * lG.nz != lE.getNs() || nt * lG.nt != lE.getNt())
	{
		throw std::invalid_argument("The lattice cannot be distributed onto the given grid.");
	}

}

LocalLatticeMemoryExtents::LocalLatticeMemoryExtents( LatticeGrid lG, LocalLatticeExtents llE, unsigned int halo_size ):
		nx(llE.nx + (lG.nx > 1 ? 2 * halo_size : 0)),
		ny(llE.ny + (lG.ny > 1 ? 2 * halo_size : 0)),
		nz(llE.nz + (lG.nz > 1 ? 2 * halo_size : 0)),
		nt(llE.nt + (lG.nt > 1 ? 2 * halo_size : 0))
{
	if(nt < halo_size) {
		throw std::invalid_argument("The lattice cannot be distributed onto the given grid.");
	}
}
