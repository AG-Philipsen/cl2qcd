/** @file
 * Declaration of the rhmc algorithm
 *
 * (c) 2013 Alessandro Sciarra <sciarra@compeng.uni-frankfurt.de>
 */

#ifndef _PHYSICS_ALGORITHMS_RHMC_
#define _PHYSICS_ALGORITHMS_RHMC_

#include "../lattices/gaugefield.hpp"
#include "../prng.hpp"
#include "../../types_hmc.h"
#include "rational_approximation.hpp"

namespace physics {
namespace algorithms {

hmc_observables perform_rhmc_step(const Rational_Approximation& approx1, const Rational_Approximation& approx2, const Rational_Approximation& approx3, const physics::lattices::Gaugefield * gf, int iter, hmc_float rnd_number, physics::PRNG& prng, const hardware::System& system);

}
}

#endif /* _PHYSICS_ALGORITHMS_RHMC_ */
