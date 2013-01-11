/** @file
 * Implementation of the solver algorithms
 *
 * (c) 2012-2013 Christopher Pinke <pinke@th.uni-frankfurt.de>
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "solver.hpp"

#include "../../logger.hpp"

static std::string create_solver_stuck_message(int iterations);

physics::algorithms::solvers::SolverStuck::SolverStuck(int iterations, std::string filename, int linenumber) : SolverException(create_solver_stuck_message(iterations), iterations, filename, linenumber) { };

static std::string create_solver_stuck_message(int iterations)
{
	std::ostringstream tmp;
	tmp << "Solver got stuck after " << iterations << " iterations";
	return tmp.str();
}
