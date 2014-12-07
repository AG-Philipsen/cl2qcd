/** @file
 * Declaration of the solver algorithms
 *
 * Copyright (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
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

#ifndef _PHYSICS_ALGORITHMS_SOLVER_
#define _PHYSICS_ALGORITHMS_SOLVER_

#include "../fermionmatrix/fermionmatrix.hpp"
#include "solvers/exceptions.hpp"
#include "solvers/cg.hpp"

namespace physics {
namespace algorithms {

  //@TODO: move everything to new file in /solver
/**
 * this namespace contains methods to solve linear systems either like \n
 * <pre>
 * A * x = b       (for x, so far only for Wilson fermions)  </pre>
 * or like
 * <pre>
 * (A + sigma) * x = b       (for x, so far only with even-odd preconditioned Staggered fields)
 * </pre>
 * @note In the latter case, the problem is solved simultaneously for an entire
 *       set of values of sigma thanks to a multi-shifted inverter
 */
namespace solvers {

/**
 * Solve the linear system A * x = b for x using the BiCGstab algorithm
 *
 *
 * \return The number of iterations performed
 *
 * \exception SolverStuck if the solver gets stuck. Contains information on performed iterations
 *
 * \exception SolverDidNotSolve if the solver did not solve (hit iteration limit).
 */
int bicgstab(const physics::lattices::Spinorfield * x, const physics::fermionmatrix::Fermionmatrix& A, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& b, const hardware::System& system, hmc_float prec);

/**
 * Solve the linear system A * x = b for x using the BiCGstab algorithm for even-odd preconditioned b and x.
 *
 *
 * \return The number of iterations performed
 *
 * \exception SolverStuck if the solver gets stuck. Contains information on performed iterations
 *
 * \exception SolverDidNotSolve if the solver did not solve (hit iteration limit).
 */
int bicgstab(const physics::lattices::Spinorfield_eo * x, const physics::fermionmatrix::Fermionmatrix_eo& A, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& b, const hardware::System& system, hmc_float prec);

}
}
}
#endif /* _PHYSICS_ALGORITHMS_SOLVER_ */
