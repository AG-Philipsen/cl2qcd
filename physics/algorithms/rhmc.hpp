/** @file
 * Declaration of the rhmc algorithm
 *
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.phys.uni-frankfurt.de>
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

#ifndef _PHYSICS_ALGORITHMS_RHMC_
#define _PHYSICS_ALGORITHMS_RHMC_

#include "../lattices/gaugefield.hpp"
#include "../prng.hpp"
#include "../../common_header_files/types_hmc.h"
#include "rational_approximation.hpp"
#include "../interfacesHandler.hpp"

namespace physics {
    namespace algorithms {

        hmc_observables perform_rhmc_step(const Rational_Approximation& approx1, const Rational_Approximation& approx2, const Rational_Approximation& approx3,
                                          const physics::lattices::Gaugefield * gf, int iter, hmc_float rnd_number, physics::PRNG& prng, const hardware::System& system,
                                          physics::InterfacesHandler& interfacesHandler);

    }
}

#endif /* _PHYSICS_ALGORITHMS_RHMC_ */
