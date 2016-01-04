/*
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

#pragma once

struct LatticeExtents
{
	LatticeExtents() : ns(4), nt(4) {};
	LatticeExtents(const int nsIn, const int ntIn): ns(nsIn), nt(ntIn) {};
	LatticeExtents(const unsigned int nsIn, const unsigned int ntIn): ns(nsIn), nt(ntIn) {};
	const unsigned int ns;
	const unsigned int nt;
};

int calculateLatticeVolume(const int nsIn, const int ntIn) noexcept;
int calculateLatticeVolume(const LatticeExtents latticeExtentsIn) noexcept;

