/** @file
 * Declaration of the heatbath algorithm
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

#ifndef _PHYSICS_ALGORITHMS_HEATBATH_
#define _PHYSICS_ALGORITHMS_HEATBATH_

#include "../prng.hpp"
#include "../lattices/gaugefield.hpp"

namespace physics {

namespace algorithms {

/**
 * Perform one heatbath step on the given gaugefield.
 *
 * Optionally includes overrelaxation.
 *
 * @param[in,out] gf The gaugefield to heatbath
 * @param[in,out] prng The PRNG to use
 * @param[in] overrelax The number of overrelaxation steps to perform. Default is none.
 */
void heatbath(physics::lattices::Gaugefield& gf, physics::PRNG& prng, int overrelax = 0);

/**
 * Perform one overrelaxation step on the given gaugefield.
 *
 * @param[in,out] gf The gaugefield to heatbath
 * @param[in,out] prng The PRNG to use
 * @param[in] steps The number of overrelaxation steps to perform. Default is 1.
 */
void overrelax(physics::lattices::Gaugefield& gf, physics::PRNG& prng, int steps = 1);

}

}

#endif /* _PHYSICS_ALGORITHMS_HEATBATH_ */
