/** @file
 * Declaration of the force algorithms
 *
 * Copyright (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
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

#ifndef _PHYSICS_ALGORITHMS_FORCES_
#define _PHYSICS_ALGORITHMS_FORCES_

#include "../lattices/gaugemomenta.hpp"
#include "../lattices/gaugefield.hpp"
#include "../lattices/spinorfield.hpp"
#include "../lattices/spinorfield_eo.hpp"
#include "../lattices/rooted_staggeredfield_eo.hpp"
#include "fermion_force.hpp"
#include "fermion_force_staggered.hpp"
#include "../interfacesHandler.hpp"

namespace physics
{
    namespace algorithms
    {

        void gauge_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Gaugefield& gf);
        void gauge_force_tlsym(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Gaugefield& gf);

        void calc_gauge_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Gaugefield& gf, physics::InterfacesHandler& interfacesHandler);

        void calc_total_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi,
                              const hardware::System& system, physics::InterfacesHandler& interfacesHandler, const physics::AdditionalParameters& additionalParameters);
        void calc_total_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi,
                              const hardware::System& system, physics::InterfacesHandler& interfacesHandler, const physics::AdditionalParameters& additionalParameters);
        void calc_total_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Gaugefield& gf, const physics::lattices::Rooted_Staggeredfield_eo& phi,
                              const hardware::System& system, physics::InterfacesHandler& interfacesHandler, const physics::AdditionalParameters& additionalParameters);

    }
}

#endif /* _PHYSICS_ALGORITHMS_FERMION_FORCES_ */
