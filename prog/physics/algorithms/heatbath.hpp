/** @file
 * Declaration of the heatbath algorithm
 */

#include "../prng.hpp"
#include "../lattices/gaugefield.hpp"

namespace physics {

namespace algorithms {

/**
 * Perform one heatbath step on the given gaugefield.
 *
 * Optionally includes overrelaxation.
 *
 * @param[in,out] gf The gaugefield to heatbath
 * @param[in,out] prng The PRNG to use
 * @param[in] overrelax The number of overrelaxation steps to perform. Default is none.
 */
void heatbath(physics::lattices::Gaugefield& gf, physics::PRNG& prng, int overrelax = 0);

/**
 * Perform one overrelaxation step on the given gaugefield.
 *
 * @param[in,out] gf The gaugefield to heatbath
 * @param[in,out] prng The PRNG to use
 * @param[in] steps The number of overrelaxation steps to perform. Default is 1.
 */
void overrelax(physics::lattices::Gaugefield& gf, physics::PRNG& prng, int steps = 1);

}

}
