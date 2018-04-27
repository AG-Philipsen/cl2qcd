/*
 * Copyright (c) 2016 Christopher Pinke
 * Copyright (c) 2016 Francesca Cuteri
 * Copyright (c) 2018 Alessandro Sciarra
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "../common_header_files/globaldefs.hpp"
#include "../hardware/size_4.hpp"  //todo: remove!
#include "latticeExtents.hpp"

typedef uint latticeIndex;

struct Index {
    Index(const latticeCoordinate x, const latticeCoordinate y, const latticeCoordinate z, const latticeCoordinate t,
          const LatticeExtents lE);
    Index(const size_4, const LatticeExtents);  // todo: remove
    Index up(const Direction dir) const;
    Index down(const Direction dir) const;
    operator latticeSize() const;
    const LatticeCoordinate x, y, z, t;
    const latticeIndex spatialIndex, globalIndex;
};

struct LinkIndex : public Index {
    LinkIndex(const Index indexIn, const Direction dirIn);
    operator latticeSize() const;
    const LinkIndex up(const Direction dirIn) const;
    const LinkIndex down(const Direction dirIn) const;
    uint get_su3_idx_ildg_format(const uint n, const uint m);  // ADD TESTS!!
    const Direction direction;
    const latticeIndex globalIndex;
};
