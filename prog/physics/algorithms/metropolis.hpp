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
#include "../lattices/gaugemomenta.hpp"
#include "../../types_hmc.h"

namespace physics {
namespace algorithms {


hmc_float calc_s_fermion(const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi, const hardware::System& system, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
hmc_float calc_s_fermion(const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);

hmc_float calc_s_fermion_mp(const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi, const hardware::System& system);
hmc_float calc_s_fermion_mp(const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system);

hmc_observables metropolis(const hmc_float rnd, const hmc_float beta, const physics::lattices::Gaugefield& gf, const physics::lattices::Gaugefield& new_u, const physics::lattices::Gaugemomenta& p, const physics::lattices::Gaugemomenta& new_p, const physics::lattices::Spinorfield& phi, const hmc_float spinor_energy_init, const hmc_float spinor_energy_mp_init, const hardware::System& system);
hmc_observables metropolis(const hmc_float rnd, const hmc_float beta, const physics::lattices::Gaugefield& gf, const physics::lattices::Gaugefield& new_u, const physics::lattices::Gaugemomenta& p, const physics::lattices::Gaugemomenta& new_p, const physics::lattices::Spinorfield_eo& phi, const hmc_float spinor_energy_init, const hmc_float spinor_energy_mp_init, const hardware::System& system);

}
}

#endif /* _PHYSICS_ALGORITHMS_METROPOLIS_ */
