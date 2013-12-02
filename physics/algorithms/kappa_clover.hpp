/** @file
 * Declaration of the kappa clover calculation
 *
 * Copyright 2012, 2013 Lars Zeidlewicz, Christopher Pinke,
 * Matthias Bach, Christian Schäfer, Stefano Lottini, Alessandro Sciarra
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

#ifndef _PHYSICS_ALGORITHMS_KAPPA_CLOVER_
#define _PHYSICS_ALGORITHMS_KAPPA_CLOVER_

#include "../prng.hpp"
#include "../lattices/gaugefield.hpp"

namespace physics {

namespace algorithms {

/**
 * Calculate kappa clover
 *
 * @param[in] gf The gaugefield to use
 * @param[in] beta
 * @return The kappa clover value
 *
 * @todo return a future
 */
hmc_float kappa_clover(physics::lattices::Gaugefield& gf, hmc_float beta);

}

}

#endif /* _PHYSICS_ALGORITHMS_KAPPA_CLOVER_ */
