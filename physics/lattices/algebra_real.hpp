/** @file
 * Declaration of some algebra operations for physics::lattices::Scalar<hmc_float>
 * together with physics::lattices::Vector<hmc_float>
 *
 * Copyright (c) 2014 Alessandro Sciarra <sciarra@th.physik.uni-frankfurt.de>
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

#ifndef _PHYSICS_LATTICES_SCALAR_REAL_
#define _PHYSICS_LATTICES_SCALAR_REAL_

#include "scalar.hpp"
#include "vector.hpp"
#include "../../common_header_files/types.h"

namespace physics {

namespace lattices {
  
//See hardware/code/real.hpp file for functions documentations
void update_zeta_cgm(const Vector<hmc_float>* out, const Vector<hmc_float>& zeta_prev, const Vector<hmc_float>& zeta_prev_prev, const Scalar<hmc_float>& sbeta_prev, const Scalar<hmc_float>& sbeta_pres, const Scalar<hmc_float>& salpha_prev, const Vector<hmc_float>& sigma, const int numeq);

void update_beta_cgm(const Vector<hmc_float>* out, const Scalar<hmc_float>& sbeta_pres, const Vector<hmc_float>& zeta_pres, const Vector<hmc_float>& zeta_prev, const int numeq);

void update_alpha_cgm(const Vector<hmc_float>* out, const Scalar<hmc_float>& salpha_pres, const Vector<hmc_float>& zeta_pres, const Vector<hmc_float>& beta_pres, const Vector<hmc_float>& zeta_prev, const Scalar<hmc_float>& sbeta_pres, const int numeq);

}

}

#endif /* _PHYSICS_LATTICES_SCALAR_REAL_ */
