/** @file
 * Implementation of the rhmc algorithm
 *
 * Copyright (c) 2013 Alessandro Sciarra <sciarra@th.phys.uni-frankfurt.de>
 *
 * This file is part of CL2QCD.
 *
 * CL2QCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * CL2QCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "rhmc.hpp"

#include "../lattices/gaugemomenta.hpp"
#include "../lattices/util.hpp"
#include "../lattices/rooted_staggeredfield_eo.hpp"
#include "../fermionmatrix/fermionmatrix_stagg.hpp"
#include "../../meta/util.hpp"
#include "integrator.hpp"
#include "metropolis.hpp"
#include "molecular_dynamics.hpp"
#include <memory>
#include "../../klepsydra/klepsydra.hpp"

// template <class SPINORFIELD> static hmc_observables perform_rhmc_step(const physics::algorithms::Rational_Approximation& approx1, const physics::algorithms::Rational_Approximation& approx2, const physics::algorithms::Rational_Approximation& approx3, const physics::lattices::Gaugefield * gf, int iter, hmc_float rnd_number, physics::PRNG& prng, const hardware::System& system);

template<class SPINORFIELD> static void init_spinorfield(const SPINORFIELD * phi, hmc_float * const spinor_energy_init, const physics::lattices::Gaugefield& gf,
        const physics::PRNG& prng, const hardware::System& system);
//Mass preconditioning not yet implemented!
//template <class SPINORFIELD> static void init_spinorfield_mp(const SPINORFIELD * phi, hmc_float * const spinor_energy_init, const SPINORFIELD * phi_mp, hmc_float * const spinor_energy_init_mp, const physics::lattices::Gaugefield& gf, const physics::PRNG& prng, const hardware::System& system);

template<class SPINORFIELD> static hmc_observables perform_rhmc_step(const physics::algorithms::Rational_Approximation& approx1,
        const physics::algorithms::Rational_Approximation& approx2, const physics::algorithms::Rational_Approximation& approx3,
        const physics::lattices::Gaugefield * const gf, const int iter, const hmc_float rnd_number, physics::PRNG& prng, const hardware::System& system,
        physics::InterfacesHandler& interfaceHandler)
{
    using namespace physics::algorithms;
    using namespace physics::lattices;

    klepsydra::Monotonic step_timer;

    const auto & params = system.get_inputparameters();

    logger.debug() << "\tRHMC:\tinit spinorfield and gaugemomentum";
    const Gaugemomenta p(system, interfaceHandler.getInterface<physics::lattices::Gaugemomenta>());
    p.gaussian(prng);

    SPINORFIELD phi(system);
    const std::auto_ptr<const SPINORFIELD> phi_mp(params.get_use_mp() ? new SPINORFIELD(system) : nullptr);
    hmc_float spinor_energy_init = 0.f;
    hmc_float spinor_energy_init_mp = 0.f;
    //Here the coefficients of phi have to be set to the rescaled ones on the base of approx1
    physics::fermionmatrix::MdagM_eo fm(system, ARG_DEF);   //with ARG_DEF, the mass is that of params
    phi.Rescale_Coefficients(approx1, fm, *gf, system, params.get_findminmax_prec(), params.get_conservative());
    if(!params.get_use_gauge_only()) {
        if(params.get_use_mp()) {
            throw Print_Error_Message("Mass preconditioning not implemented for staggered fermions!", __FILE__, __LINE__);
            //init_spinorfield_mp(&phi, &spinor_energy_init, phi_mp.get(), &spinor_energy_init_mp, *gf, prng, system);
        } else {
            init_spinorfield(&phi, &spinor_energy_init, *gf, prng, system);
        }
    }

    logger.debug() << "\tRHMC:\tupdate gaugefield and gaugemomentum";
    const GaugefieldParametersImplementation gaugefieldParameters { &params };
    const Gaugefield new_u(system, &gaugefieldParameters, prng, false);
    const Gaugemomenta new_p(system, interfaceHandler.getInterface<physics::lattices::Gaugemomenta>());
    // copy u->u' p->p' for the integrator
    copyData(&new_u, *gf);
    copyData(&new_p, p);

    //here, clmem_phi is inverted several times and stored in clmem_phi_inv
    logger.debug() << "\tRHMC:\tcall integrator";
    //Before MD the coefficients of phi have to be set to the rescaled ones on the base of approx2
    phi.Rescale_Coefficients(approx2, fm, *gf, system, params.get_findminmax_prec(), params.get_conservative());
    if(params.get_use_mp()) {
        throw Print_Error_Message("Mass preconditioning not implemented for staggered fermions!", __FILE__, __LINE__);
        //integrator(&new_p, &new_u, phi, *phi_mp.get(), system);
    } else {
        integrator(&new_p, &new_u, phi, system, interfaceHandler);
    }

    //metropolis step: afterwards, the updated config is again in gaugefield and p
    logger.debug() << "\tRHMC [MET]:\tperform Metropolis step: ";
    //Before Metropolis test the coeff. of phi have to be set to the rescaled ones on the base of approx3
    phi.Rescale_Coefficients(approx3, fm, *gf, system, params.get_findminmax_prec(), params.get_conservative());
    //this call calculates also the HMC-Observables
    const hmc_observables obs = metropolis(rnd_number, params.get_beta(), *gf, new_u, p, new_p, phi, spinor_energy_init, phi_mp.get(), spinor_energy_init_mp,
            system, interfaceHandler);

    if(obs.accept == 1) {
        // perform the change nonprimed->primed !
        copyData(gf, new_u);
        logger.info() << "\tRHMC [MET]:\tnew configuration accepted";
    } else {
        logger.info() << "\tRHMC [MET]:\tnew configuration rejected";
    }
    logger.info() << "\tRHMC:\tfinished trajectory " << iter;
    logger.info() << "\tRHMC:\tstep duration (ms): " << step_timer.getTime() / 1e3f;

    return obs;
}

hmc_observables physics::algorithms::perform_rhmc_step(const physics::algorithms::Rational_Approximation& approx1,
        const physics::algorithms::Rational_Approximation& approx2, const physics::algorithms::Rational_Approximation& approx3,
        const physics::lattices::Gaugefield * const gf, const int iter, const hmc_float rnd_number, physics::PRNG& prng, const hardware::System& system,
        physics::InterfacesHandler& interfaceHandler)
{
    using namespace physics::lattices;

    const auto & params = system.get_inputparameters();
    if(params.get_use_eo()) {
        return ::perform_rhmc_step<Rooted_Staggeredfield_eo>(approx1, approx2, approx3, gf, iter, rnd_number, prng, system, interfaceHandler);
    } else {
        throw Print_Error_Message("RHMC algorithm not implemented for non even-odd preconditioned fields!", __FILE__, __LINE__);
        //return ::perform_rhmc_step<Rooted_Staggeredfield>(approx1, approx2, approx3, gf, iter, rnd_number, prng, system);
    }
}

template<class SPINORFIELD> static void init_spinorfield(const SPINORFIELD * phi, hmc_float * const spinor_energy_init, const physics::lattices::Gaugefield& gf,
        const physics::PRNG& prng, const hardware::System& system)
{
    using namespace physics::algorithms;

    const auto & params = system.get_inputparameters();

    const SPINORFIELD initial(system);

    //init/update spinorfield phi
    initial.set_gaussian(prng);
    //calc init energy for spinorfield
    *spinor_energy_init = squarenorm(initial);
    //update spinorfield
    md_update_spinorfield(phi, gf, initial, system, params.get_mass());
}

//Mass preconditioning not yet implemented!
/*
 template <class SPINORFIELD> static void init_spinorfield_mp(const SPINORFIELD * phi, hmc_float * const spinor_energy_init, const SPINORFIELD * phi_mp, hmc_float * const spinor_energy_init_mp, const physics::lattices::Gaugefield& gf, const physics::PRNG& prng, const hardware::System& system)
 {
 using namespace physics::algorithms;

 const auto & params = system.get_inputparameters();

 const SPINORFIELD initial(system);

 //init/update spinorfield phi
 initial.gaussian(prng);
 //calc init energy for spinorfield
 *spinor_energy_init = squarenorm(initial);
 //update spinorfield with heavy mass: det(kappa_mp, mu_mp)
 md_update_spinorfield(phi, gf, initial, system, params.get_kappa_mp(), meta::get_mubar_mp(params));
 initial.gaussian(prng);
 //calc init energy for mass-prec spinorfield (this is the same as for the spinorfield above)
 *spinor_energy_init_mp = squarenorm(initial);
 //update detratio spinorfield: det(kappa, mu) / det(kappa_mp, mu_mp)
 md_update_spinorfield_mp(phi_mp, gf, initial, system);
 }
 */
