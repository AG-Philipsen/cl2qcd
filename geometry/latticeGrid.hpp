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

#pragma once

#include "latticeExtents.hpp"
#include "../common_header_files/globaldefs.h"
#include <ostream>

struct FourUnsignedInt
{
	FourUnsignedInt(const unsigned int x, const unsigned int y, const unsigned int z, const unsigned int t);
	const unsigned int x, y, z, t;
};
std::ostream& operator<<(std::ostream& out, const FourUnsignedInt & integers);


struct LatticeGrid : public FourUnsignedInt
{
	LatticeGrid( const unsigned int numberOfDevices);
};

struct LatticeGridIndex : public FourUnsignedInt
{
	LatticeGridIndex(const unsigned int x, const unsigned int y, const unsigned int z, const unsigned int t, const LatticeGrid);
};

struct LocalLatticeExtents : public FourUnsignedInt
{
	LocalLatticeExtents(const LatticeGrid lG, const LatticeExtents lE);
};

struct LocalLatticeMemoryExtents : public FourUnsignedInt
{
	LocalLatticeMemoryExtents(const LatticeGrid lG, const LocalLatticeExtents llE, unsigned int halo_size);
};
