/** @file
 * Implementation of the physics::lattices::Rooted_Staggeredfield_eo class
 *
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
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
	
	Set_a0(aux.Get_a0());
	Set_a(aux.Get_a());
	Set_b(aux.Get_b());
}


