/**
 * Copyright 2016 Christopher Pinke
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
# pragma once

#include "latticeGrid.hpp"

//todo: introduce similar class which is for Non-Link types (different hypervolume)
struct TemporalParallelizationHandler
{
	TemporalParallelizationHandler(const LatticeGridIndex, const LocalLatticeExtents, const size_t sizeOfElement, const size_t haloSize);
	latticeIndex getMainPartIndex_destination() const;
	latticeIndex getFirstHaloIndex_destination() const;
	latticeIndex getSecondHaloIndex_destination() const;
	size_t getMainPartSize() const;
	size_t getHaloPartSize() const ;
	size_t getMainPartSizeInBytes() const;
	size_t getHaloPartSizeInBytes() const ;
	latticeIndex getMainPartIndex_source() const;
	latticeIndex getFirstHaloPartIndex_source() const;
	latticeIndex getSecondHaloPartIndex_source() const;
	virtual latticeSize hyperVolume() const = 0;
	const LocalLatticeExtents lLE;
	const LatticeGridIndex lGI;
	const latticeSize localExtentInSlowestDirection;
	const latticeSize haloSize;
	const size_t sizeInBytesPerElement;
	const Direction slowestDirection;
};

struct TemporalParallelizationHandlerLink : public TemporalParallelizationHandler
{
	TemporalParallelizationHandlerLink(const LatticeGridIndex, const LocalLatticeExtents, const size_t sizeOfElement, const size_t haloSize);
	latticeSize hyperVolume() const override;
};

struct TemporalParallelizationHandlerNonLink : public TemporalParallelizationHandler
{
	TemporalParallelizationHandlerNonLink(const LatticeGridIndex, const LocalLatticeExtents, const size_t sizeOfElement, const size_t haloSize);
	latticeSize hyperVolume() const override;
};