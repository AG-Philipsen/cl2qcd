/** @file
 * Implementation of the physics::lattices::Rooted_Staggeredfield_eo class
 * 
 * (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
 */

#include <algorithm>    // std::max
#include "rooted_staggeredfield_eo.hpp"


// #include "../../meta/util.hpp"
// #include <cassert>
// #include "../../hardware/code/spinors_staggered.hpp"
// #include "../../hardware/code/fermions_staggered.hpp"
// #include "../../meta/type_ops.hpp"
// #include "../../hardware/buffers/halo_update.hpp"
// //For hardware::code::get_eoprec_spinorfieldsize()
// #include "../../hardware/code/spinors.hpp"


physics::lattices::Rooted_Staggeredfield_eo::Rooted_Staggeredfield_eo(const hardware::System& system)
	: Staggeredfield_eo(system), physics::algorithms::Rational_Coefficients(std::max(system.get_inputparameters().get_metro_approx_ord(), system.get_inputparameters().get_md_approx_ord()))
{
}

physics::lattices::Rooted_Staggeredfield_eo::Rooted_Staggeredfield_eo(const physics::algorithms::Rational_Approximation& approx, const hardware::System& system)
	: Staggeredfield_eo(system), physics::algorithms::Rational_Coefficients(approx.Get_order(), approx.Get_a0(), approx.Get_a(), approx.Get_b()) 
{
}

void physics::lattices::Rooted_Staggeredfield_eo::Rescale_Coefficients(const physics::algorithms::Rational_Approximation& approx, const physics::fermionmatrix::Fermionmatrix_stagg_eo& A, const physics::lattices::Gaugefield& gf, const hardware::System& system, hmc_float prec, bool conservative)
{
	physics::algorithms::Rational_Coefficients aux = approx.Rescale_Coefficients(A, gf, system, prec, conservative);
	
	Set_coeff(aux.Get_a0(), aux.Get_a(), aux.Get_b());
}


