/** @file
 * Declaration of the inversion algorithms
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

#ifndef _PHYSICS_ALGORITHMS_INVERSION_
#define _PHYSICS_ALGORITHMS_INVERSION_

#include "../lattices/gaugefield.hpp"
#include "../lattices/spinorfield.hpp"

namespace physics {

    namespace algorithms {

        /**
         * Perform the inversion and store result to solution_buffer
         *
         * @param[out] result Spinorfield in which to store the inversion result
         * @param[in] gaugefield Gaugefield on which to base the inversion
         * @param[in] sources Spinorfields from which to start the inversion
         */
        void perform_inversion(const std::vector<physics::lattices::Spinorfield*> * result, const physics::lattices::Gaugefield* gaugefield,
                               const std::vector<physics::lattices::Spinorfield*>& sources, const hardware::System& system,
                               physics::InterfacesHandler& interfacesHandler);

    }

}

#endif /* _PHYSICS_ALGORITHMS_INVERSION_ */
