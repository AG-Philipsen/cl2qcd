/** @file
 * Declaration of the kappa clover calculation
 */

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
