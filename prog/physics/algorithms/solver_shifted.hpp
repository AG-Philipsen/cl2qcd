/** @file
 * Declaration of the solver algorithms
 *
 * (c) 2013 Alessandro Sciarra <sciarra@th.uni-frankfurt.de>
 */

#ifndef _PHYSICS_ALGORITHMS_SOLVER_SHIFTED_
#define _PHYSICS_ALGORITHMS_SOLVER_SHIFTED_

#include "../fermionmatrix/fermionmatrix_stagg.hpp"

namespace physics {
namespace algorithms {
namespace solvers {

/**
 * Solve the linear system (A + sigma) * x = b for x and a whole set of values of sigma,
 * using the CG-M algorithm (multi-shifted inverter) for even-odd preconditioned b and x.
 * For further details, see B.Jegerlehner arXiv:hep-lat/9612014v1, even if here there are some missing information
 * and also some wrong indeces. Then compare with "Accurate conjugate gradient methods for families of
 * shifted systems" that is a work by Jasper van den Eshof and Gerard L. G. Sleijpen of December 2003:
 * here the notation is completely different, but one can build the complete CG-M algorithm. This was done
 * and summarized in the file CG-M.pdf (see Feature #482)
 * 
 * \return The number of iterations performed
 *
 * \exception SolverStuck if the solver gets stuck. Contains information on performed iterations
 *
 * \exception SolverDidNotSolve if the solver did not solve (hit iteration limit).
 */
int cg_m(const std::vector<physics::lattices::Staggeredfield_eo *> x, const std::vector<hmc_float> sigma, const physics::fermionmatrix::Fermionmatrix_stagg_eo& A, const physics::lattices::Gaugefield& gf, const physics::lattices::Staggeredfield_eo& b, const hardware::System& system, hmc_float prec);

}
}
}
#endif /* _PHYSICS_ALGORITHMS_SOLVER_SHIFTED_ */
