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

#ifndef _PHYSICS_ALGORITHMS_SOLVERS_
#define _PHYSICS_ALGORITHMS_SOLVERS_

#include "../../fermionmatrix/fermionmatrix.hpp"
#include "exceptions.hpp"
#include "../../lattices/util.hpp"
#include "../../lattices/scalar_complex.hpp"

#include "cg.hpp"
#include "bicgstab.hpp"

namespace physics {
namespace algorithms {

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


}
}
}

#endif
