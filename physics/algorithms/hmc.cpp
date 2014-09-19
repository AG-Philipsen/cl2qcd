/** @file
 * Declaration of the hmc algorithm
 *
 * Copyright (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 * Copyright (c) 2012-2013 Christopher Pinke <pinke@compeng.uni-frankfurt.de>
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

#include "hmc.hpp"

#include "../lattices/gaugemomenta.hpp"
#include "../lattices/util.hpp"
#include "../lattices/spinorfield.hpp"
#include "../lattices/spinorfield_eo.hpp"
#include "../../meta/util.hpp"
#include "integrator.hpp"
#include "metropolis.hpp"
#include "molecular_dynamics.hpp"
#include <memory>
#include "../../klepsydra/klepsydra.hpp"

template <class SPINORFIELD> static hmc_observables perform_hmc_step(const physics::lattices::Gaugefield * gf, int iter, hmc_float rnd_number, physics::PRNG& prng, const hardware::System& system);

template <class SPINORFIELD> static void init_spinorfield(const SPINORFIELD * phi, hmc_float * const spinor_energy_init, const physics::lattices::Gaugefield& gf, const physics::PRNG& prng, const hardware::System& system);
template <class SPINORFIELD> static void init_spinorfield_mp(const SPINORFIELD * phi, hmc_float * const spinor_energy_init, const SPINORFIELD * phi_mp, hmc_float * const spinor_energy_init_mp, const physics::lattices::Gaugefield& gf, const physics::PRNG& prng, const hardware::System& system);


template <class SPINORFIELD> static hmc_observables perform_hmc_step(const physics::lattices::Gaugefield * const gf, const int iter, const hmc_float rnd_number, physics::PRNG& prng, const hardware::System& system)
{
	using namespace physics::algorithms;
	using namespace physics::lattices;

	klepsydra::Monotonic step_timer;

	const auto & params = system.get_inputparameters();

	logger.trace() << "\tHMC:\tinit spinorfield and gaugemomentum" ;
	const Gaugemomenta p(system);
	p.gaussian(prng);

	const SPINORFIELD phi(system);
	const std::auto_ptr<const SPINORFIELD> phi_mp(params.get_use_mp() ? new SPINORFIELD(system) : nullptr);
	hmc_float spinor_energy_init = 0.f;
	hmc_float spinor_energy_init_mp = 0.f;
	if(!params.get_use_gauge_only()) {
		if(params.get_use_mp()) {
			init_spinorfield_mp(&phi, &spinor_energy_init, phi_mp.get(), &spinor_energy_init_mp, *gf, prng, system);
		} else {
			init_spinorfield(&phi, &spinor_energy_init, *gf, prng, system);
		}
	}

	logger.trace() << "\tHMC:\tupdate gaugefield and gaugemomentum" ;
	const Gaugefield new_u(system, prng, false);
	const Gaugemomenta new_p(system);
	// copy u->u' p->p' for the integrator
	copyData(&new_u, *gf);
	copyData(&new_p, p);

	//here, clmem_phi is inverted several times and stored in clmem_phi_inv
	logger.trace() << "\tHMC:\tcall integrator" ;
	if(params.get_use_mp()) {
		integrator(&new_p, &new_u, phi, *phi_mp.get(), system);
	} else {
		integrator(&new_p, &new_u, phi, system);
	}

	//metropolis step: afterwards, the updated config is again in gaugefield and p
	logger.trace() << "\tHMC [MET]:\tperform Metropolis step: ";
	//this call calculates also the HMC-Observables
	const hmc_observables obs = metropolis(rnd_number, params.get_beta(), *gf, new_u, p, new_p, phi, spinor_energy_init, phi_mp.get(), spinor_energy_init_mp, system);

	if(obs.accept == 1) {
		// perform the change nonprimed->primed !
		copyData(gf, new_u);
		logger.info() << "\tHMC [MET]:\tnew configuration accepted" ;
	} else {
		logger.info() << "\tHMC [MET]:\tnew configuration rejected" ;
	}
	logger.info() << "\tHMC:\tfinished trajectory " << iter ;
	logger.info() << "\tHMC:\tstep duration (ms): " << step_timer.getTime() / 1e3f;

	return obs;
}

hmc_observables physics::algorithms::perform_hmc_step(const physics::lattices::Gaugefield * const gf, const int iter, const hmc_float rnd_number, physics::PRNG& prng, const hardware::System& system)
{
	using namespace physics::lattices;

	const auto & params = system.get_inputparameters();
	if(params.get_use_eo()) {
		return ::perform_hmc_step<Spinorfield_eo>(gf, iter, rnd_number, prng, system);
	} else {
		return ::perform_hmc_step<Spinorfield>(gf, iter, rnd_number, prng, system);
	}
}

template <class SPINORFIELD> static void init_spinorfield(const SPINORFIELD * phi, hmc_float * const spinor_energy_init, const physics::lattices::Gaugefield& gf, const physics::PRNG& prng, const hardware::System& system)
{
	using namespace physics::algorithms;

	const auto & params = system.get_inputparameters();

	const SPINORFIELD initial(system);

	//init/update spinorfield phi
	initial.gaussian(prng);
	//calc init energy for spinorfield
	*spinor_energy_init = squarenorm(initial);
	//update spinorfield: det(kappa, mu)
	md_update_spinorfield(phi, gf, initial, system, params.get_kappa(), meta::get_mubar(params));
}

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
