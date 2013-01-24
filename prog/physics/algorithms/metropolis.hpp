/** @file
 * Declaration of the metropolis algorithm
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */


#ifndef _PHYSICS_ALGORITHMS_METROPOLIS_
#define _PHYSICS_ALGORITHMS_METROPOLIS_

#include "../lattices/gaugefield.hpp"
#include "../lattices/spinorfield.hpp"
#include "../lattices/spinorfield_eo.hpp"
#include "../../types_hmc.h"

namespace physics {
namespace algorithms {


hmc_float calc_s_fermion(const physics::lattices::Spinorfield * phi_inv, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi, const hardware::System& system, hmc_float kappa, hmc_float mubar);
hmc_float calc_s_fermion(const physics::lattices::Spinorfield_eo * phi_inv, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, hmc_float kappa, hmc_float mubar);

hmc_float calc_s_fermion_mp(const physics::lattices::Spinorfield * phi_inv, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi, const hardware::System& system, hmc_float kappa, hmc_float mubar);
hmc_float calc_s_fermion_mp(const physics::lattices::Spinorfield_eo * phi_inv, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, hmc_float kappa, hmc_float mubar);

}
}

#endif /* _PHYSICS_ALGORITHMS_METROPOLIS_ */
