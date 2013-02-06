/** @file
 * Declaration of the force algorithms
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#ifndef _PHYSICS_ALGORITHMS_FORCES_
#define _PHYSICS_ALGORITHMS_FFORCES_

#include "../lattices/gaugemomenta.hpp"
#include "../lattices/gaugefield.hpp"
#include "../lattices/spinorfield.hpp"
#include "../lattices/spinorfield_eo.hpp"
#include "fermion_force.hpp"

namespace physics {
namespace algorithms {

void gauge_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Gaugefield& gf);
void gauge_force_tlsym(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Gaugefield& gf);

void calc_gauge_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Gaugefield& gf, const hardware::System& system);

void calc_total_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi, const hardware::System& system, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
void calc_total_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);


}
}

#endif /* _PHYSICS_ALGORITHMS_FERMION_FORCES_ */
