/**
 * Copyright 2015 Christopher Pinke
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

int calculateLatticeVolume(const int nsIn, const int ntIn) noexcept
{
	return 	nsIn * nsIn * nsIn * ntIn;
}

int calculateLatticeVolume(const LatticeExtents latticeExtentsIn) noexcept
{
	return 	calculateLatticeVolume(latticeExtentsIn.ns, latticeExtentsIn.nt);
}

int calculateSpatialLatticeVolume(const int nsIn) noexcept
{
	return nsIn * nsIn * nsIn;
}

int calculateSpatialLatticeVolume(const LatticeExtents latticeExtentsIn) noexcept
{
	return calculateSpatialLatticeVolume(latticeExtentsIn.ns);
}



