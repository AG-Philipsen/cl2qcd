/** @file
 * Declaration of some basic operations for physics::lattices::Scalar<hmc_complex>
 *
 * Copyright (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
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

#include "scalar.hpp"
#include "../../common_header_files/types.h"

namespace physics {

namespace lattices {

void add(const Scalar<hmc_complex>* dest, const Scalar<hmc_complex>& left, const Scalar<hmc_complex>& right);
void subtract(const Scalar<hmc_complex>* dest, const Scalar<hmc_complex>& minuend, const Scalar<hmc_complex>& subtrahend);
void multiply(const Scalar<hmc_complex>* dest, const Scalar<hmc_complex>& left, const Scalar<hmc_complex>& right);
void divide(const Scalar<hmc_complex>* dest, const Scalar<hmc_complex>& numerator, const Scalar<hmc_complex>& denominator);

void convert(const Scalar<hmc_complex>* dest, const Scalar<hmc_float>& src);

}

}
