/** @file
 * Declaration of the molecular dynamics algorithms
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

#ifndef _PHYSICS_ALGORITHMS_MOLECULAR_DYNAMICS_
#define _PHYSICS_ALGORITHMS_MOLECULAR_DYNAMICS_

#include "../lattices/gaugefield.hpp"
#include "../lattices/spinorfield.hpp"
#include "../lattices/spinorfield_eo.hpp"
#include "../lattices/rooted_staggeredfield_eo.hpp"
#include "../lattices/gaugemomenta.hpp"
#include "rational_approximation.hpp"
#include "../interfacesHandler.hpp"

namespace physics {

    namespace algorithms {

        void md_update_gaugefield(const physics::lattices::Gaugefield * gf, const physics::lattices::Gaugemomenta&, hmc_float eps);

        void md_update_spinorfield(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf,
                                   const physics::lattices::Spinorfield& orig, const hardware::System& system,
                                   physics::InterfacesHandler & interfacesHandler, const physics::AdditionalParameters& additionalParameters);
        void md_update_spinorfield(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf,
                                   const physics::lattices::Spinorfield_eo& orig, const hardware::System& system,
                                   physics::InterfacesHandler & interfacesHandler, const physics::AdditionalParameters& additionalParameters);
        void md_update_spinorfield(const physics::lattices::Rooted_Staggeredfield_eo * out, const physics::lattices::Gaugefield& gf,
                                   const physics::lattices::Rooted_Staggeredfield_eo& orig, const hardware::System& system,
                                   physics::InterfacesHandler & interfacesHandler, const physics::AdditionalParameters& additionalParameters);
        void md_update_spinorfield_mp(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf,
                                      const physics::lattices::Spinorfield& orig, const hardware::System& system, physics::InterfacesHandler& interfacesHandler);
        void md_update_spinorfield_mp(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf,
                                      const physics::lattices::Spinorfield_eo& orig, const hardware::System& system, physics::InterfacesHandler& interfacesHandler);

        void md_update_gaugemomentum(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf,
                                     const physics::lattices::Spinorfield& phi, const hardware::System& system,
                                     physics::InterfacesHandler& interfaceHandler, const physics::AdditionalParameters& additionalParameters);
        void md_update_gaugemomentum(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf,
                                     const physics::lattices::Spinorfield_eo& phi, const hardware::System& system,
                                     physics::InterfacesHandler& interfaceHandler, const physics::AdditionalParameters& additionalParameters);
        void md_update_gaugemomentum(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf,
                                     const physics::lattices::Rooted_Staggeredfield_eo& phi, const hardware::System& system,
                                     physics::InterfacesHandler& interfaceHandler, const physics::AdditionalParameters& additionalParameters);

        void md_update_gaugemomentum_gauge(const physics::lattices::Gaugemomenta * const gm, hmc_float eps, const physics::lattices::Gaugefield& gf,
                                           const hardware::System& system, physics::InterfacesHandler& interfaceHandler);
        void md_update_gaugemomentum_fermion(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf,
                                             const physics::lattices::Spinorfield& phi, const hardware::System& system,
                                             physics::InterfacesHandler& interfaceHandler, const physics::AdditionalParameters& additionalParameters);
        void md_update_gaugemomentum_fermion(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf,
                                             const physics::lattices::Spinorfield_eo& phi, const hardware::System& system,
                                             physics::InterfacesHandler& interfaceHandler, const physics::AdditionalParameters& additionalParameters);
        void md_update_gaugemomentum_fermion(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf,
                                             const physics::lattices::Rooted_Staggeredfield_eo& phi, const hardware::System& system,
                                             physics::InterfacesHandler& interfaceHandler, const physics::AdditionalParameters& additionalParameters);
        void md_update_gaugemomentum_detratio(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf,
                                              const physics::lattices::Spinorfield& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler);
        void md_update_gaugemomentum_detratio(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf,
                                              const physics::lattices::Spinorfield_eo& phi_mp, const hardware::System& system, physics::InterfacesHandler& interfaceHandler);
    }

}

#endif /* _PHYSICS_ALGORITHMS_MOLECULAR_DYNAMICS_ */
