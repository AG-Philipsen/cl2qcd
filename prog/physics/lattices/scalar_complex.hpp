/** @file
 * Declaration of some basic operations for physics::lattices::Scalar<hmc_complex>
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#include "scalar.hpp"
#include "../../types.h"

namespace physics {

namespace lattices {

void add(const Scalar<hmc_complex>* dest, const Scalar<hmc_complex>& left, const Scalar<hmc_complex>& right);
void subtract(const Scalar<hmc_complex>* dest, const Scalar<hmc_complex>& minuend, const Scalar<hmc_complex>& subtrahend);
void multiply(const Scalar<hmc_complex>* dest, const Scalar<hmc_complex>& left, const Scalar<hmc_complex>& right);
void divide(const Scalar<hmc_complex>* dest, const Scalar<hmc_complex>& numerator, const Scalar<hmc_complex>& denominator);

void convert(const Scalar<hmc_complex>* dest, const Scalar<hmc_float>& src);

}

}
