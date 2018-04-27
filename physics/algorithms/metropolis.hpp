/** @file
 * Declaration of the metropolis algorithm
 *
 * Copyright (c) 2013,2015,2016,2018 Alessandro Sciarra
 * Copyright (c) 2013 Matthias Bach
 * Copyright (c) 2015 Christopher Czaban
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

#ifndef _PHYSICS_ALGORITHMS_METROPOLIS_
#define _PHYSICS_ALGORITHMS_METROPOLIS_

#include "../../common_header_files/types_hmc.hpp"
#include "../interfacesHandler.hpp"
#include "../lattices/gaugefield.hpp"
#include "../lattices/gaugemomenta.hpp"
#include "../lattices/rooted_staggeredfield_eo.hpp"
#include "../lattices/spinorfield.hpp"
#include "../lattices/spinorfield_eo.hpp"

namespace physics {
    namespace algorithms {

        hmc_float calc_s_fermion(const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi,
                                 const hardware::System& system, physics::InterfacesHandler& interfacesHandler,
                                 const physics::AdditionalParameters& additionalParameters);
        hmc_float calc_s_fermion(const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi,
                                 const hardware::System& system, physics::InterfacesHandler& interfacesHandler,
                                 const physics::AdditionalParameters& additionalParameters);
        hmc_float calc_s_fermion(const physics::lattices::Gaugefield& gf,
                                 const physics::lattices::Rooted_Staggeredfield_eo& phi, const hardware::System& system,
                                 physics::InterfacesHandler& interfacesHandler,
                                 const physics::AdditionalParameters& additionalParameters);

        hmc_float calc_s_fermion_mp(const physics::lattices::Gaugefield& gf,
                                    const physics::lattices::Spinorfield& phi_mp, const hardware::System& system,
                                    physics::InterfacesHandler& interfacesHandler);
        hmc_float calc_s_fermion_mp(const physics::lattices::Gaugefield& gf,
                                    const physics::lattices::Spinorfield_eo& phi_mp, const hardware::System& system,
                                    physics::InterfacesHandler& interfacesHandler);
        hmc_float
        calc_s_fermion_mp(const physics::lattices::Gaugefield& gf,
                          const physics::lattices::Rooted_Staggeredfield_eo& phi_mp, const hardware::System& system,
                          physics::InterfacesHandler& interfacesHandler);  // function so far NOT IMPLEMENTED!!

        hmc_observables metropolis(const hmc_float rnd, const hmc_float beta, const physics::lattices::Gaugefield& gf,
                                   const physics::lattices::Gaugefield& new_u, const physics::lattices::Gaugemomenta& p,
                                   const physics::lattices::Gaugemomenta& new_p,
                                   const physics::lattices::Spinorfield& phi, const hmc_float spinor_energy_init,
                                   const physics::lattices::Spinorfield* const phi_mp,
                                   const hmc_float spinor_energy_mp_init, const hardware::System& system,
                                   physics::InterfacesHandler& interfacesHandler);
        hmc_observables metropolis(const hmc_float rnd, const hmc_float beta, const physics::lattices::Gaugefield& gf,
                                   const physics::lattices::Gaugefield& new_u, const physics::lattices::Gaugemomenta& p,
                                   const physics::lattices::Gaugemomenta& new_p,
                                   const physics::lattices::Spinorfield_eo& phi, const hmc_float spinor_energy_init,
                                   const physics::lattices::Spinorfield_eo* const phi_mp,
                                   const hmc_float spinor_energy_mp_init, const hardware::System& system,
                                   physics::InterfacesHandler& interfacesHandler);
        hmc_observables
        metropolis(const hmc_float rnd, const hmc_float beta, const physics::lattices::Gaugefield& gf,
                   const physics::lattices::Gaugefield& new_u, const physics::lattices::Gaugemomenta& p,
                   const physics::lattices::Gaugemomenta& new_p, const physics::lattices::Rooted_Staggeredfield_eo& phi,
                   const hmc_float spinor_energy_init, const physics::lattices::Rooted_Staggeredfield_eo* const phi_mp,
                   const hmc_float spinor_energy_mp_init, const hardware::System& system,
                   physics::InterfacesHandler& interfacesHandler);  // mass preconditioning is so far NOT IMPLEMENTED!!

    }  // namespace algorithms
}  // namespace physics

#endif /* _PHYSICS_ALGORITHMS_METROPOLIS_ */
