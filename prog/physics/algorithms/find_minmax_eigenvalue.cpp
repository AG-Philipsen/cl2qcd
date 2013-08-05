/** @file
 * Implementation of the algorithm to find min and max eigenvalues of an operator
 * 
 * (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
 */

#include "find_minmax_eigenvalue.hpp"
#include "../lattices/staggeredfield_eo.hpp"
#include "../lattices/util.hpp"
#include "../lattices/scalar_complex.hpp"


hmc_float find_max_eigenvalue(const physics::fermionmatrix::Fermionmatrix_stagg_eo& A, const physics::lattices::Gaugefield& gf, const hardware::System& system, hmc_float prec)
{
	using namespace physics::lattices;
	using namespace physics::algorithms;
	
	hmc_float max;
	const Scalar<hmc_complex> norm(system);
	const Scalar<hmc_complex> norm_prev(system);
	Staggeredfield_eo v(system);
	pseudo_randomize<Staggeredfield_eo, su3vec>(&v, 123);
	
}



