/** @file
 * Declaration of the flavour doublets algorithms
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

#ifndef _PHYSICS_ALGORITHMS_FLAVOUR_DOUBLET_
#define _PHYSICS_ALGORITHMS_FLAVOUR_DOUBLET_

#include "../lattices/gaugefield.hpp"
#include "../lattices/spinorfield.hpp"

namespace physics {

namespace algorithms {

/**
 * Calculate flavour multiplet correlators from private solution_buffer and store them to a file
 * @param[in] gaugefield The gaugefield
 * @param[in] result The Spinorfields
 * @param[in] corr_fn filename
 */
void flavour_doublet_correlators(const std::vector<physics::lattices::Spinorfield*>& result, const std::vector<physics::lattices::Spinorfield*>& sources, std::string corr_fn, const hardware::System& system);

std::vector<hmc_float> calculate_correlator(const std::string& type, const std::vector<physics::lattices::Spinorfield*>& corr, const std::vector<physics::lattices::Spinorfield*>& sources, const hardware::System& system);
}

}

#endif /* _PHYSICS_ALGORITHMS_FLAVOUR_DOUBLET_ */
