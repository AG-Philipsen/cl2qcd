/** @file
 * Declaration of the integrator algorithms
 *
 * See hep-lat/0505020 for more infos on integrators.
 * Currently, two Integrators are implemented: Leapfrog and 2MN.
 * Leapfrog-Integration is the simplest second-order integration, thus the erros made are not
 *  small with regard to 2MN-integration. However, it only requires one force-calculation per step.
 * 2MN (2order Minimized Norm) needs two force-calculations per step, but has a much lower error then leapfrog. It also
 *  features one tunable parameter, lambda.
 * If one has the Hamiltonian:
 *  H = p^2/2 + S(q) = T + V
 * then one leapfrog step of length eps (= tau/number_of_steps) is
 *  exp(eps(T + V) ) = exp(eps/2 T) exp( eps V ) exp( eps/2 T)
 * and one 2MN-step is
 *  exp(eps(T + V) ) = exp(lambda*eps T) exp( eps/2 V ) exp( (1 - 2lamdba) *eps T) exp( eps/2 V ) exp( lamdba*eps T)
 * One can also use more than one timescale for the different pars to improve speed (e.g. hep-lat/0209037 and hep-lat/0506011v2).
 * Up to now, 2 timescales are implemented (hardcoded). The Hamiltonian is split up
 *  H = T + V_gauge + V_fermion = S_1 + S_2
 * In case of the Leapfrog this means then
 *  exp(eps(S_1) ) = exp(eps/2 T(V_gauge)) exp( eps V ) exp( eps/2 T (V_gauge) )
 *  exp(eps(S_2) ) = exp(eps T(V_fermion))
 *  -> exp(eps(T + V) ) = exp(eps/2 S_1) [ exp( eps/m S_2 ) ]^m exp( eps/2 S_1)
 * In case of the 2MN this means then
 *  exp(eps(S_1) ) = exp(lambda*eps T(V_gauge) ) exp( eps/2 V ) exp( (1 - 2lamdba) *eps T(V_gauge) ) exp( eps/2 V ) exp( lamdba*eps T(V_gauge) )
 *  exp(eps(S_2) ) = exp(eps T(V_fermion))
 *  -> exp(eps(T + V) ) = exp(lambda*eps/2 S_1) [ exp( eps/m S_2 ) ]^m exp( (1 - 2lamdba)*eps/2 S_1) [ exp( eps/m S_2 ) ]^m exp( lamdba*eps/2 S_1)
 * In the program:
 *  exp(eps T) == md_update_gaugemomentum(eps)
 *  exp(eps V) == md_update_gaugefield(eps)
 * NOTE:  Performing number_of_steps integrationsteps leads to "whole" steps for the momentum
 *          ( exp(eps/2 T)exp(eps/2 T) = exp(eps T) )
 *
 * Copyright (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 * Copyright (c) 2012-2013 Christopher Pinke <pinke@compeng.uni-frankfurt.de>
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

#ifndef _PHYSICS_ALGORITHMS_INTEGRATOR_
#define _PHYSICS_ALGORITHMS_INTEGRATOR_

#include "../lattices/gaugemomenta.hpp"
#include "../lattices/gaugefield.hpp"
#include "../lattices/spinorfield.hpp"
#include "../lattices/spinorfield_eo.hpp"
#include "../lattices/rooted_staggeredfield_eo.hpp"
#include "../interfacesHandler.hpp"

namespace physics {

    namespace algorithms {

        void integrator(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
                        const physics::lattices::Spinorfield& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler);
        void integrator(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
                        const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler);
        void integrator(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
                        const physics::lattices::Rooted_Staggeredfield_eo& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler);
        void integrator(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
                        const physics::lattices::Spinorfield& phi, const physics::lattices::Spinorfield& phi_mp, const hardware::System& system,
                        physics::InterfacesHandler& interfaceHandler);
        void integrator(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
                        const physics::lattices::Spinorfield_eo& phi, const physics::lattices::Spinorfield_eo& phi_mp, const hardware::System& system,
                        physics::InterfacesHandler& interfaceHandler);

        void leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
                      const physics::lattices::Spinorfield& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler);
        void leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
                      const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler);
        void leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
                      const physics::lattices::Rooted_Staggeredfield_eo& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler);
        void leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
                      const physics::lattices::Spinorfield& phi, const physics::lattices::Spinorfield& phi_mp, const hardware::System& system,
                      physics::InterfacesHandler& interfaceHandler);
        void leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
                      const physics::lattices::Spinorfield_eo& phi, const physics::lattices::Spinorfield_eo& phi_mp, const hardware::System& system,
                      physics::InterfacesHandler& interfaceHandler);

        void twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const physics::lattices::Spinorfield& phi,
                   const hardware::System& system, physics::InterfacesHandler& interfaceHandler);
        void twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
                   const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler);
        void twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
                   const physics::lattices::Rooted_Staggeredfield_eo& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler);
        void twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const physics::lattices::Spinorfield& phi,
                   const physics::lattices::Spinorfield& phi_mp, const hardware::System& system, physics::InterfacesHandler& interfaceHandler);
        void twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
                   const physics::lattices::Spinorfield_eo& phi, const physics::lattices::Spinorfield_eo& phi_mp, const hardware::System& system,
                   physics::InterfacesHandler& interfaceHandler);

    }

}

#endif /* _PHYSICS_ALGORITHMS_INTEGRATOR_ */
