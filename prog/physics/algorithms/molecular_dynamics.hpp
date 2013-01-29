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

void gauge_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Gaugefield& gf);
void gauge_force_tlsym(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Gaugefield& gf);

void fermion_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Spinorfield& Y, const physics::lattices::Spinorfield& X, const physics::lattices::Gaugefield& gf, hmc_float kappa = ARG_DEF);
void fermion_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Spinorfield_eo& Y, const physics::lattices::Spinorfield_eo& X, int evenodd, const physics::lattices::Gaugefield& gf, hmc_float kappa = ARG_DEF);

void calc_fermion_forces(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi, const hardware::System& system, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
void calc_fermion_forces(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);

void calc_gauge_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Gaugefield& gf, const hardware::System& system);

void calc_total_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi, const hardware::System& system, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
void calc_total_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);

void md_update_gaugemomentum(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi, const hardware::System& system, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
void md_update_gaugemomentum(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);

}

}

#endif /* _PHYSICS_ALGORITHMS_MOLECULAR_DYNAMICS_ */
