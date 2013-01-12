/** @file
 * Declaration of the solver algorithms
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#ifndef _PHYSICS_ALGORITHMS_SOLVER_
#define _PHYSICS_ALGORITHMS_SOLVER_

#include "../fermionmatrix/fermionmatrix.hpp"
#include "../../exceptions.h"

namespace physics {
namespace algorithms {

/**
 * this namespace contains methods to solve linear systems like
 * A * x = b
 * for x
 */
namespace solvers {

/**
 * Base exception used by solvers ot indicate solving failures.
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

/**
 * Thrown by solvers to indicate being stuck.
 */
class SolverStuck : public SolverException {
public:
	SolverStuck(int iterations, std::string filename, int linenumber);
};

/**
 * Thrown by solvers to indicate that it could not solve
 */
class SolverDidNotSolve : public SolverException {
public:
	SolverDidNotSolve(int iterations, std::string filename, int linenumber) : SolverException("Solver did not solve.", iterations, filename, linenumber) { };
};

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

/**
 * Solve the linear system A * x = b for x using the CG algorithm
 *
 *
 * \return The number of iterations performed
 *
 * \exception SolverStuck if the solver gets stuck. Contains information on performed iterations
 *
 * \exception SolverDidNotSolve if the solver did not solve (hit iteration limit).
 */
int cg(const physics::lattices::Spinorfield * x, const physics::fermionmatrix::Fermionmatrix& A, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& b, const hardware::System& system, hmc_float prec);

}
}
}
#endif /* _PHYSICS_ALGORITHMS_SOLVER_ */
