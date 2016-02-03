/** @file
 * Declaration of the fermion_force functions
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

#ifndef _PHYSICS_ALGORITHMS_FERMION_FORCE_STAGGERED_
#define _PHYSICS_ALGORITHMS_FERMION_FORCE_STAGGERED_

#include "../lattices/gaugemomenta.hpp"
#include "../lattices/gaugefield.hpp"
#include "../lattices/rooted_staggeredfield_eo.hpp"
#include "rational_approximation.hpp"
#include "../interfacesHandler.hpp"

namespace physics
{
    namespace algorithms
    {

        //These methods really calculate the total fermion force and they add it to the Gaugemomenta field
        void calc_fermion_forces(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf,
                                 const physics::lattices::Rooted_Staggeredfield_eo& phi, const hardware::System& system,
                                 physics::InterfacesHandler& interfacesHandler, const physics::AdditionalParameters& additionalParameters);

        //Here, in the following functions, there is the detailed force calculation
        void calc_fermion_force(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf,
                                const physics::lattices::Rooted_Staggeredfield_eo& phi, const hardware::System& system,
                                physics::InterfacesHandler& interfacesHandler, const physics::AdditionalParameters& additionalParameters);

        //These methods interfaces only the lower level of the code (Molecular_Dynamics class) with the upper one,
        //namely they just call the function that enqueues the kernel
        void fermion_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Staggeredfield_eo& A,
                           const physics::lattices::Staggeredfield_eo& B, const physics::lattices::Gaugefield& gf, int evenodd);

    }
}

#endif /* _PHYSICS_ALGORITHMS_FERMION_FORCE_STAGGERED_ */
