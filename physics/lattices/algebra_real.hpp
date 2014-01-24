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
/*
void update_zeta_cgm(const Scalar<hmc_float>* dest, const Scalar<hmc_float>& left, const Scalar<hmc_complex>& right);

const hardware::buffers::Plain<hmc_float> * zeta_prev, const hardware::buffers::Plain<hmc_float> * zeta_prev_prev, const hardware::buffers::Plain<hmc_float> * sbeta_prev, const hardware::buffers::Plain<hmc_float> * sbeta_pres, const hardware::buffers::Plain<hmc_float> * salpha_prev, const hardware::buffers::Plain<hmc_float> * sigma, const int numeq, const hardware::buffers::Plain<hmc_float> * out



void subtract(const Scalar<hmc_complex>* dest, const Scalar<hmc_complex>& minuend, const Scalar<hmc_complex>& subtrahend);
void multiply(const Scalar<hmc_complex>* dest, const Scalar<hmc_complex>& left, const Scalar<hmc_complex>& right);
void divide(const Scalar<hmc_complex>* dest, const Scalar<hmc_complex>& numerator, const Scalar<hmc_complex>& denominator);

void convert(const Scalar<hmc_complex>* dest, const Scalar<hmc_float>& src);
*/
}

}

#endif /* _PHYSICS_LATTICES_SCALAR_REAL_ */
