/** @file
 * Declaration of some algebra operations for physics::lattices::Scalar<hmc_float>
 * together with physics::lattices::Vector<hmc_float>
 *
 * Copyright (c) 2014,2018 Alessandro Sciarra
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef _PHYSICS_LATTICES_SCALAR_REAL_
#define _PHYSICS_LATTICES_SCALAR_REAL_

#include "scalar.hpp"
#include "vector.hpp"
#include "../../common_header_files/types.hpp"
#include "../../hardware/system.hpp"

namespace physics {

namespace lattices {

///@todo These two functions should be replaced by more general ones in the Vector class!!!
void access_real_vector_element(const Scalar<hmc_float>* out, const Vector<hmc_float>& in, const int index);
void access_real_vector_element(const Vector<hmc_float>* out, const Scalar<hmc_float>& in, const int index);

void add(const Scalar<hmc_float>* dest, const Scalar<hmc_float>& left, const Scalar<hmc_float>& right);
void subtract(const Scalar<hmc_float>* dest, const Scalar<hmc_float>& minuend, const Scalar<hmc_float>& subtrahend);
void multiply(const Scalar<hmc_float>* dest, const Scalar<hmc_float>& left, const Scalar<hmc_float>& right);
void divide(const Scalar<hmc_float>* dest, const Scalar<hmc_float>& numerator, const Scalar<hmc_float>& denominator);

//See hardware/code/real.hpp file for functions documentations
void update_zeta_cgm(const Vector<hmc_float>* out, const Vector<hmc_float>& zeta_prev, const Vector<hmc_float>& zeta_prev_prev, const Scalar<hmc_float>& sbeta_prev, const Scalar<hmc_float>& sbeta_pres, const Scalar<hmc_float>& salpha_prev, const Vector<hmc_float>& sigma, const int numeq);

void update_beta_cgm(const Vector<hmc_float>* out, const Scalar<hmc_float>& sbeta_pres, const Vector<hmc_float>& zeta_pres, const Vector<hmc_float>& zeta_prev, const int numeq);

void update_alpha_cgm(const Vector<hmc_float>* out, const Scalar<hmc_float>& salpha_pres, const Vector<hmc_float>& zeta_pres, const Vector<hmc_float>& beta_pres, const Vector<hmc_float>& zeta_prev, const Scalar<hmc_float>& sbeta_pres, const int numeq);

size_t get_flops_update_cgm(const std::string quantity, const int Neqs, const hardware::System&);

}

}

#endif /* _PHYSICS_LATTICES_SCALAR_REAL_ */
