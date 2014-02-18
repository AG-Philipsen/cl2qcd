/** @file
 * Declaration of the staggered chiral condensate calculation.
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

#ifndef _PHYSICS_ALGORITHMS_CHIRAL_CONDENSATE_STAGG_
#define _PHYSICS_ALGORITHMS_CHIRAL_CONDENSATE_STAGG_

#include "../lattices/gaugefield.hpp"

namespace physics {

namespace algorithms {

/**
 * Calculate chiral condesate with noise estimators (only volume source is accepted so far)
 * with the content and the number of the sources as specified in Inputparameters
 * @param[in] gf The actual configuration used in the fermionmatrix
 * @param[in] prng The actual random number generator
 * @param[in] system The system to operate on
 */
hmc_complex chiral_condensate_staggered(const physics::lattices::Gaugefield& gf, physics::PRNG& prng, const hardware::System& system);

}

}

#endif /* _PHYSICS_ALGORITHMS_CHIRAL_CONDENSATE_STAGG_ */
