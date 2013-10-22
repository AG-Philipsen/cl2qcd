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
#include "../lattices/rooted_staggeredfield_eo.hpp"
#include "../lattices/gaugemomenta.hpp"
#include "rational_approximation.hpp"

namespace physics {

namespace algorithms {

void md_update_gaugefield(const physics::lattices::Gaugefield * gf, const physics::lattices::Gaugemomenta& , hmc_float eps);
void md_update_spinorfield(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& orig, const hardware::System& system, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
void md_update_spinorfield(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& orig, const hardware::System& system, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
void md_update_spinorfield(const physics::lattices::Staggeredfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Staggeredfield_eo& orig, const physics::algorithms::Rational_Coefficients& coeff, const hardware::System& system, hmc_float mass);
void md_update_spinorfield_mp(const physics::lattices::Spinorfield * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& orig, const hardware::System& system, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
void md_update_spinorfield_mp(const physics::lattices::Spinorfield_eo * out, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& orig, const hardware::System& system, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);

void md_update_gaugemomentum(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi, const hardware::System& system, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
void md_update_gaugemomentum(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
void md_update_gaugemomentum(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf, const physics::lattices::Rooted_Staggeredfield_eo& phi, const hardware::System& system, hmc_float mass, hmc_float mubar = ARG_DEF); //mubar is so far NOT IMPLEMENTED!!

void md_update_gaugemomentum_gauge(const physics::lattices::Gaugemomenta * const gm, hmc_float eps, const physics::lattices::Gaugefield& gf, const hardware::System& system);
void md_update_gaugemomentum_fermion(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi, const hardware::System& system, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
void md_update_gaugemomentum_fermion(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
void md_update_gaugemomentum_fermion(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf, const physics::lattices::Rooted_Staggeredfield_eo& phi, const hardware::System& system, hmc_float mass, hmc_float mubar = ARG_DEF); //mubar is so far NOT IMPLEMENTED!!
void md_update_gaugemomentum_detratio(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi, const hardware::System& system);
void md_update_gaugemomentum_detratio(const physics::lattices::Gaugemomenta * const inout, hmc_float eps, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi_mp, const hardware::System& system);
}

}

#endif /* _PHYSICS_ALGORITHMS_MOLECULAR_DYNAMICS_ */
