/** @file
 * Declaration of the algorithm to find min and max eigenvalues of an operator
 * 
 * (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
 */

#ifndef _PHYSICS_ALGORITHMS_FINDMINMAX_
#define _PHYSICS_ALGORITHMS_FINDMINMAX_

#include "../lattices/gaugefield.hpp"
#include "../fermionmatrix/fermionmatrix_stagg.hpp"

namespace physics {

namespace algorithms {

/**
 * This function returns the maximum eigenvalue of the operator A up to some specified precision.
 *  @param A The operator whose maximum eigenvalue has to be found
 *  @param gf The gaugefield on which A depends
 *  @param system The system it is operating on
 *  @param prec The precision up to which the maximum eigenvalue is found
 *
 *  @note The algorithm here implemented is the so called "Power Method". Hence to be sure
 *        that the algorithm works, the eigenvalues must be real. If this is not the case,
 *        the algorithm could give meaningless results or not converge. Actually, even if
 *        the eigenvalues are real, if there is not only ONE dominant eigenvalue, then 
 *        the algorithm won't converge (for example if there are two opposite biggest eigenvalues).
 *        However, to be honest, in the RHMC with staggered fermions, one needs to find the
 *        maximum and the minimum eigenvalues of the matrix MdagM that is hermitian and then with
 *        real positive eigenvalue. For this reason, this function will calculate the eigenvalue
 *        if and only if the operator A is hermitian.
 *  
 */
hmc_float find_max_eigenvalue(const physics::fermionmatrix::Fermionmatrix_stagg_eo& A, const physics::lattices::Gaugefield& gf, const hardware::System& system, hmc_float prec);

}

}

#endif /* _PHYSICS_ALGORITHMS_FINDMINMAX_ */
