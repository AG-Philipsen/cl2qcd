/** @file
 * Implementation of the rhmc algorithm
 *
 * Copyright (c) 2013, 2017 Alessandro Sciarra <sciarra@th.phys.uni-frankfurt.de>
 * Copyright (c) 2017 Francesca Cuteri <cuteri@th.physik.uni-frankfurt.de>
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

template<class SPINORFIELD> static void init_spinorfield(const SPINORFIELD * phi, hmc_float * const spinor_energy_init, const physics::lattices::Gaugefield& gf,
        const physics::PRNG& prng, const hardware::System& system, physics::InterfacesHandler& interfacesHandler);

template<class SPINORFIELD> static hmc_observables perform_rhmc_step(const physics::algorithms::Rational_Approximation& approx1,
        const physics::algorithms::Rational_Approximation& approx2, const physics::algorithms::Rational_Approximation& approx3,
        const physics::lattices::Gaugefield * const gf, const int iter, const hmc_float rnd_number, physics::PRNG& prng, const hardware::System& system,
        physics::InterfacesHandler& interfacesHandler)
{
    using namespace physics::algorithms;
    using namespace physics::lattices;

    klepsydra::Monotonic step_timer;

    const physics::algorithms::RhmcParametersInterface & parametersInterface = interfacesHandler.getRhmcParametersInterface();
    const physics::AdditionalParameters& additionalParameters = interfacesHandler.getAdditionalParameters<SPINORFIELD>();

    logger.debug() << "\tRHMC:\tinit spinorfield and gaugemomentum";
    const Gaugemomenta p(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());
    p.gaussian(prng);

    SPINORFIELD phi(system, interfacesHandler.getInterface<SPINORFIELD>());
    const std::auto_ptr<const SPINORFIELD> phi_mp(parametersInterface.getUseMp() ? new SPINORFIELD(system, interfacesHandler.getInterface<SPINORFIELD>()) : nullptr);
    hmc_float spinor_energy_init = 0.f;
    hmc_float spinor_energy_init_mp = 0.f;
    //Here the coefficients of phi have to be set to the rescaled ones on the base of approx1
    physics::fermionmatrix::MdagM_eo fm(system, interfacesHandler.getInterface<physics::fermionmatrix::MdagM_eo>());
    hmc_float maxEigenvalue;
    hmc_float minEigenvalue;
    find_maxmin_eigenvalue(maxEigenvalue, minEigenvalue, fm, *gf, system, interfacesHandler, parametersInterface.getFindMinMaxPrec(), additionalParameters);
    if(parametersInterface.getConservative())
        maxEigenvalue *= 1.05;
    hmc_float conditionNumber=maxEigenvalue/minEigenvalue;
	phi.Rescale_Coefficients(approx1, minEigenvalue, maxEigenvalue);

    if(!parametersInterface.getUseGaugeOnly()) {
        if(parametersInterface.getUseMp()) {
            throw Print_Error_Message("Mass preconditioning not implemented for staggered fermions!", __FILE__, __LINE__);
            //init_spinorfield_mp(&phi, &spinor_energy_init, phi_mp.get(), &spinor_energy_init_mp, *gf, prng, system);
        } else {
            init_spinorfield(&phi, &spinor_energy_init, *gf, prng, system, interfacesHandler);
        }
    }

    logger.debug() << "\tRHMC:\tupdate gaugefield and gaugemomentum";
    const Gaugefield new_u(system, &interfacesHandler.getInterface<physics::lattices::Gaugefield>(), prng, false);
    const Gaugemomenta new_p(system, interfacesHandler.getInterface<physics::lattices::Gaugemomenta>());
    // copy u->u' p->p' for the integrator
    copyData(&new_u, *gf);
    copyData(&new_p, p);

    //here, clmem_phi is inverted several times and stored in clmem_phi_inv
    logger.debug() << "\tRHMC:\tcall integrator";
    //Before MD the coefficients of phi have to be set to the rescaled ones on the base of approx2
    phi.Rescale_Coefficients(approx2, minEigenvalue, maxEigenvalue);
    if(parametersInterface.getUseMp()) {
        throw Print_Error_Message("Mass preconditioning not implemented for staggered fermions!", __FILE__, __LINE__);
        //integrator(&new_p, &new_u, phi, *phi_mp.get(), system);
    } else {
        integrator(&new_p, &new_u, phi, system, interfacesHandler);
    }

    //metropolis step: afterwards, the updated config is again in gaugefield and p
    logger.debug() << "\tRHMC [MET]:\tperform Metropolis step: ";
    //Before Metropolis test the coeff. of phi have to be set to the rescaled ones on the base of approx3
    find_maxmin_eigenvalue(maxEigenvalue, minEigenvalue, fm, new_u, system, interfacesHandler, parametersInterface.getFindMinMaxPrec(), additionalParameters);
	if(parametersInterface.getConservative())
	    maxEigenvalue *= 1.05;
	phi.Rescale_Coefficients(approx3, minEigenvalue, maxEigenvalue);
    //this call calculates also the HMC-Observables
    const hmc_observables obs = metropolis(rnd_number, parametersInterface.getBeta(), *gf, new_u, p, new_p, phi, spinor_energy_init,
                                           phi_mp.get(), spinor_energy_init_mp, system, interfacesHandler);

    if(obs.accept == 1) {
        // perform the change nonprimed->primed !
        copyData(gf, new_u);
        logger.info() << "\tRHMC [MET]:\tnew configuration accepted";
    } else {
        logger.info() << "\tRHMC [MET]:\tnew configuration rejected";
    }
    logger.info() << "\tRHMC:\tfinished trajectory " << iter;
    logger.info() << "\tRHMC:\tstep duration (ms): " << step_timer.getTime() / 1e3f;

	//The optimal number of pseudofermion is calculated according to http://arxiv.org/abs/hep-lat/0608015v1
	std::ostringstream reportOnConditionNumber;
	reportOnConditionNumber << "\tRHMC:\tcondition number for trajectory " << iter << ": lambda_max/lambda_min = " << conditionNumber;
	unsigned optimalNumberPseudofermions = std::floor(0.5*std::log(conditionNumber));
	if(optimalNumberPseudofermions <= 1)
	    reportOnConditionNumber << "  =>  multiple pseudofermions method not advantageous!";
	else
	    reportOnConditionNumber << "  =>  multiple pseudofermions method suggested  =>  optimal number of pseudofermion: n_opt = " << optimalNumberPseudofermions;
	logger.info() << reportOnConditionNumber.str();

    return obs;
}

hmc_observables physics::algorithms::perform_rhmc_step(const physics::algorithms::Rational_Approximation& approx1,
        const physics::algorithms::Rational_Approximation& approx2, const physics::algorithms::Rational_Approximation& approx3,
        const physics::lattices::Gaugefield * const gf, const int iter, const hmc_float rnd_number, physics::PRNG& prng, const hardware::System& system,
        physics::InterfacesHandler& interfaceHandler)
{
    using namespace physics::lattices;

    const physics::algorithms::RhmcParametersInterface & parametersInterface = interfaceHandler.getRhmcParametersInterface();
    if(parametersInterface.getUseEo()) {
        return ::perform_rhmc_step<Rooted_Staggeredfield_eo>(approx1, approx2, approx3, gf, iter, rnd_number, prng, system, interfaceHandler);
    } else {
        throw Print_Error_Message("RHMC algorithm not implemented for non even-odd preconditioned fields!", __FILE__, __LINE__);
        //return ::perform_rhmc_step<Rooted_Staggeredfield>(approx1, approx2, approx3, gf, iter, rnd_number, prng, system);
    }
}

template<class SPINORFIELD> static void init_spinorfield(const SPINORFIELD * phi, hmc_float * const spinor_energy_init, const physics::lattices::Gaugefield& gf,
        const physics::PRNG& prng, const hardware::System& system, physics::InterfacesHandler& interfacesHandler)
{
    using namespace physics::algorithms;

    const physics::AdditionalParameters& additionalParameters = interfacesHandler.getAdditionalParameters<SPINORFIELD>();
    const SPINORFIELD initial(system, interfacesHandler.getInterface<SPINORFIELD>());

    //init/update spinorfield phi
    initial[0]->set_gaussian(prng);
    //calc init energy for spinorfield
    *spinor_energy_init = squarenorm(*initial[0]);
    //update spinorfield
    md_update_spinorfield(phi, gf, initial, system, interfacesHandler, additionalParameters);
}

