/** @file
 * Declaration of the metropolis algorithm
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

#ifndef _PHYSICS_ALGORITHMS_FERMION_FORCE_
#define _PHYSICS_ALGORITHMS_FERMION_FORCE_

#include "../lattices/gaugemomenta.hpp"
#include "../lattices/gaugefield.hpp"
#include "../lattices/spinorfield.hpp"
#include "../lattices/spinorfield_eo.hpp"
#include "../lattices/staggeredfield_eo.hpp"
#include "rational_approximation.hpp"

namespace physics {
namespace algorithms {

//These methods really calculate the total fermion force and they add it to the Gaugemomenta field
void calc_fermion_forces(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi, const hardware::System& system, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
void calc_fermion_forces(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
void calc_fermion_forces(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf, const physics::lattices::Staggeredfield_eo& phi, const physics::algorithms::Rational_Coefficients& coeff, const hardware::System& system, hmc_float mass);

void calc_detratio_forces(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi, const hardware::System& system);
void calc_detratio_forces(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system);

//Here, in the following functions, there is the detailed force calculation (these functions
//are called from those above, that are actually unified by a template in the .cpp file,
//except for the staggered case)
void calc_fermion_force(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi, const hardware::System& system, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
void calc_fermion_force(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
void calc_fermion_force(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf, const physics::lattices::Staggeredfield_eo& phi, const physics::algorithms::Rational_Coefficients& coeff, const hardware::System& system, hmc_float mass);

void calc_fermion_force_detratio(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield& phi_mp, const hardware::System& system);
void calc_fermion_force_detratio(const physics::lattices::Gaugemomenta * force, const physics::lattices::Gaugefield& gf, const physics::lattices::Spinorfield_eo& phi_mp, const hardware::System& system);

//These methods interfaces only the lower level of the code (Molecular_Dynamics class) with the upper one,
//namely the just call the function that enqueues the kernel
void fermion_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Spinorfield& Y, const physics::lattices::Spinorfield& X, const physics::lattices::Gaugefield& gf, hmc_float kappa = ARG_DEF);
void fermion_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Spinorfield_eo& Y, const physics::lattices::Spinorfield_eo& X, int evenodd, const physics::lattices::Gaugefield& gf, hmc_float kappa = ARG_DEF);
void fermion_force(const physics::lattices::Gaugemomenta * gm, const physics::lattices::Staggeredfield_eo& A, const physics::lattices::Staggeredfield_eo& X, int evenodd);

}
}

#endif /* _PHYSICS_ALGORITHMS_FERMION_FORCE_ */
