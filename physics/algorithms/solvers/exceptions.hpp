/** @file
 * Declaration of solver exceptions
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

#ifndef _PHYSICS_ALGORITHMS_SOLVER_EXCEPTIONS_
#define _PHYSICS_ALGORITHMS_SOLVER_EXCEPTIONS_

#include "../../executables/exceptions.h"

namespace physics {
namespace algorithms {
namespace solvers {

/**
 * Base exception used by solvers to indicate solving failures.
 */
class SolverException : public Print_Error_Message {
public:
	/**
	 * How many iterations the solver performed until the iteration occured.
	 */
	int get_iterations() const noexcept {
		return iterations;
	};
protected:
	/**
	 * Create a solver exception with a printable message and the iteration at which it occured.
	 *
	 * Protected to not allow creation of generic solver exceptions.
	 */
	SolverException(std::string message, int iterations, std::string filename, int linenumber) : Print_Error_Message(message, filename, linenumber), iterations(iterations) { };
private:
	const int iterations;
};

class SolverStuck : public SolverException {
public:
	SolverStuck(int iterations, std::string filename, int linenumber);
};

class SolverDidNotSolve : public SolverException {
public:
	SolverDidNotSolve(int iterations, std::string filename, int linenumber) : SolverException("Solver did not solve.", iterations, filename, linenumber) { };
};

}
}
}

#endif
