/*
 * Copyright 2016 Francesca Cuteri, Christopher Pinke
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
#include "index.hpp"

struct LatticeGrid : public LatticeExtents
{
	LatticeGrid(const uint numberOfDevices, const LatticeExtents);
};

struct LatticeGridIndex : public Index
{
	LatticeGridIndex(const latticeCoordinate x, const latticeCoordinate y, const latticeCoordinate z, const latticeCoordinate t, const LatticeGrid);
};

struct LocalLatticeExtents: public LatticeExtents
{
	LocalLatticeExtents(const LatticeExtents, const LatticeGrid);
};

struct LocalLatticeMemoryExtents: public LatticeExtents
{
	LocalLatticeMemoryExtents(const LatticeGrid, const LocalLatticeExtents, uint);
};

