/** @file
 * Declaration of the bicgstab algorithms
 *
 * Copyright (c) 2012,2013 Matthias Bach
 * Copyright (c) 2012,2014 Christopher Pinke
 * Copyright (c) 2015,2016,2018 Alessandro Sciarra
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

#ifndef _PHYSICS_ALGORITHMS_SOLVERS_BICGSTAB_
#define _PHYSICS_ALGORITHMS_SOLVERS_BICGSTAB_

#include "../../interfacesHandler.hpp"
#include "solvers.hpp"

namespace physics {
    namespace algorithms {
        namespace solvers {

            /**
             * Solve the linear system A * x = b for x using the BiCGstab algorithm
             *
             * \return The number of iterations performed
             * \exception SolverStuck if the solver gets stuck. Contains information on performed iterations
             * \exception SolverDidNotSolve if the solver did not solve (hit iteration limit).
             */
            int bicgstab(const physics::lattices::Spinorfield* x, const physics::fermionmatrix::Fermionmatrix& A,
                         const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& b,
                         const hardware::System& system, physics::InterfacesHandler& interfacesHandler, hmc_float prec,
                         const physics::AdditionalParameters& additionalParameters);

            /**
             * Solve the linear system A * x = b for x using the BiCGstab algorithm for even-odd preconditioned b and x.
             *
             * \return The number of iterations performed
             * \exception SolverStuck if the solver gets stuck. Contains information on performed iterations
             * \exception SolverDidNotSolve if the solver did not solve (hit iteration limit).
             */
            int bicgstab(const physics::lattices::Spinorfield_eo* x, const physics::fermionmatrix::Fermionmatrix_eo& A,
                         const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& b,
                         const hardware::System& system, physics::InterfacesHandler& interfacesHandler, hmc_float prec,
                         const physics::AdditionalParameters& additionalParameters);

        }  // namespace solvers
    }      // namespace algorithms
}  // namespace physics
#endif
