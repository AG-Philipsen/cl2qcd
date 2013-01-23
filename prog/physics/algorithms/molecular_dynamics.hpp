/** @file
 * Declaration of the molecular dynamics algorithms
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#ifndef _PHYSICS_ALGORITHMS_MOLECULAR_DYNAMICS_
#define _PHYSICS_ALGORITHMS_MOLECULAR_DYNAMICS_

#include "../lattices/gaugefield.hpp"
#include "../lattices/spinorfield.hpp"
#include "../lattices/spinorfield_eo.hpp"
#include "../lattices/gaugemomenta.hpp"

namespace physics {

namespace algorithms {

void md_update_gaugemomenta(const physics::lattices::Gaugemomenta * dest, const physics::lattices::Gaugemomenta& src, hmc_float eps);
void md_update_gaugefield(const physics::lattices::Gaugefield * gf, const physics::lattices::Gaugemomenta& , hmc_float eps);
void gauge_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Gaugefield& gf);
void gauge_force_tlsym(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Gaugefield& gf);
void fermion_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Spinorfield& Y, const physics::lattices::Spinorfield& X, const physics::lattices::Gaugefield& gf, hmc_float kappa = ARG_DEF);
void fermion_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Spinorfield_eo& Y, const physics::lattices::Spinorfield_eo& X, int evenodd, const physics::lattices::Gaugefield& gf, hmc_float kappa = ARG_DEF);

}

}

#endif /* _PHYSICS_ALGORITHMS_MOLECULAR_DYNAMICS_ */
