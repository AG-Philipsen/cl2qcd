/** @file
 * Declaration of the kappa clover calculation
 */

#ifndef _PHYSICS_ALGORITHMS_KAPPA_CLOVER_
#define _PHYSICS_ALGORITHMS_KAPPA_CLOVER_

#include "../prng.hpp"
#include "../lattices/gaugefield.hpp"

namespace physics {

namespace algorithms {

/**
 * Calculate kappa clover
 *
 * @param[in] gf The gaugefield to use
 * @param[in] beta
 * @return The kappa clover value
 *
 * @todo return a future
 */
hmc_float kappa_clover(physics::lattices::Gaugefield& gf, hmc_float beta);

}

}

#endif /* _PHYSICS_ALGORITHMS_KAPPA_CLOVER_ */
