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

#include "parallelization.hpp"

TemporalParallelizationHandler::TemporalParallelizationHandler(const LatticeGridIndex lGI, const LocalLatticeExtents lLE, const size_t sizeOfElement, const size_t haloSize) :
	lLE(lLE),
	lGI(lGI),
	localExtentInSlowestDirection(lLE.tExtent),
	haloSize(haloSize),
	sizeInBytesPerElement(sizeOfElement),
	slowestDirection(TDIR)
{}

TemporalParallelizationHandlerLink::TemporalParallelizationHandlerLink(const LatticeGridIndex lGI, const LocalLatticeExtents lLE, const size_t sizeOfElement, const size_t haloSize):
	TemporalParallelizationHandler(lGI, lLE, sizeOfElement, haloSize)
{}

latticeSize TemporalParallelizationHandlerLink::hyperVolume() const
{
	return lLE.xExtent*lGI.x.extent * lLE.yExtent*lGI.y.extent * lLE.zExtent*lGI.z.extent * NDIM;
}

TemporalParallelizationHandlerNonLink::TemporalParallelizationHandlerNonLink(const LatticeGridIndex lGI, const LocalLatticeExtents lLE, const size_t sizeOfElement, const size_t haloSize):
	TemporalParallelizationHandler(lGI, lLE, sizeOfElement, haloSize)
{}

latticeSize TemporalParallelizationHandlerNonLink::hyperVolume() const
{
	return lLE.xExtent*lGI.x.extent * lLE.yExtent*lGI.y.extent * lLE.zExtent*lGI.z.extent;
}

latticeIndex TemporalParallelizationHandler::getMainPartIndex_destination() const
{
	return 0;
}

latticeIndex TemporalParallelizationHandler::getFirstHaloIndex_destination() const
{
	return hyperVolume()*localExtentInSlowestDirection;
}

latticeIndex TemporalParallelizationHandler::getSecondHaloIndex_destination() const
{
	return hyperVolume()*(localExtentInSlowestDirection + haloSize);
}

size_t TemporalParallelizationHandler::getMainPartSize() const
{
	return hyperVolume()*localExtentInSlowestDirection;
}

size_t TemporalParallelizationHandler::getHaloPartSize() const
{
	return hyperVolume()*haloSize;
}

size_t TemporalParallelizationHandler::getMainPartSizeInBytes() const
{
	return hyperVolume()*localExtentInSlowestDirection * sizeInBytesPerElement;
}

size_t TemporalParallelizationHandler::getHaloPartSizeInBytes() const
{
	return hyperVolume()*haloSize * sizeInBytesPerElement;
}

latticeIndex TemporalParallelizationHandler::getMainPartIndex_source() const
{
	return lGI * hyperVolume()*localExtentInSlowestDirection;
}

latticeIndex TemporalParallelizationHandler::getFirstHaloPartIndex_source() const
{
	return lGI.up(slowestDirection) * hyperVolume()*localExtentInSlowestDirection;
}

latticeIndex TemporalParallelizationHandler::getSecondHaloPartIndex_source() const
{
	return hyperVolume() * ((lGI * localExtentInSlowestDirection - haloSize + lLE.tExtent * lGI.t.extent) % (lLE.tExtent * lGI.t.extent));
}


