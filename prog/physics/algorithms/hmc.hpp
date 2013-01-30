/** @file
 * Declaration of the hmc algorithm
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 * (c) 2012-2013 Christopher Pinke <pinke@compeng.uni-frankfurt.de>
 */

#ifndef _PHYSICS_ALGORITHMS_HMC_
#define _PHYSICS_ALGORITHMS_HMC_

#include "../lattices/gaugefield.hpp"
#include "../prng.hpp"

namespace physics {
namespace algorithms {

hmc_observables perform_hmc_step(const physics::lattices::Gaugefield * gf, int iter, hmc_float rnd_number, physics::PRNG& prng, const hardware::System& system);

}
}

#endif /* _PHYSICS_ALGORITHMS_HMC_ */
