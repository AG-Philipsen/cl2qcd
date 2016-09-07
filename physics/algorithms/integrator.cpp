/** @file
 * Implementation of the integrator algorithms
 *
 * Copyright (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 * Copyright (c) 2012-2013 Christopher Pinke <pinke@compeng.uni-frankfurt.de>
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

#include "integrator.hpp"
#include "molecular_dynamics.hpp"
#include "../../meta/util.hpp"

template<class SPINORFIELD> static void integrator(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const SPINORFIELD& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler);
template<class SPINORFIELD> static void integrator(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const SPINORFIELD& phi, const SPINORFIELD& phi_mp, const hardware::System& system, physics::InterfacesHandler& interfaceHandler);

template<class SPINORFIELD> static void leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const SPINORFIELD& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler);
template<class SPINORFIELD> static void leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const SPINORFIELD& phi, const SPINORFIELD& phi_mp, const hardware::System& system, physics::InterfacesHandler& interfaceHandler);
template<class SPINORFIELD> static void leapfrog_1ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const SPINORFIELD& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler);
template<class SPINORFIELD> static void leapfrog_2ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const SPINORFIELD& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler);
template<class SPINORFIELD> static void leapfrog_3ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const SPINORFIELD& phi, const SPINORFIELD& phi_mp, const hardware::System& system, physics::InterfacesHandler& interfaceHandler);

template<class SPINORFIELD> static void twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const SPINORFIELD& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler);
template<class SPINORFIELD> static void twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const SPINORFIELD& phi, const SPINORFIELD& phi_mp, const hardware::System& system, physics::InterfacesHandler& interfaceHandler);
template<class SPINORFIELD> static void twomn_1ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const SPINORFIELD& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler);
template<class SPINORFIELD> static void twomn_2ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const SPINORFIELD& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler);
template<class SPINORFIELD> static void twomn_3ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const SPINORFIELD& phi, const SPINORFIELD& phi_mp, const hardware::System& system, physics::InterfacesHandler& interfaceHandler);

static void check_integrator_params(physics::InterfacesHandler& interfaceHandler);

void physics::algorithms::leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const physics::lattices::Spinorfield& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler)
{
    ::leapfrog(gm, gf, phi, system, interfaceHandler);
}
void physics::algorithms::leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler)
{
    ::leapfrog(gm, gf, phi, system, interfaceHandler);
}
void physics::algorithms::leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const physics::lattices::Rooted_Staggeredfield_eo& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler)
{
    ::leapfrog(gm, gf, phi, system, interfaceHandler);
}
void physics::algorithms::leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const physics::lattices::Spinorfield& phi, const physics::lattices::Spinorfield& phi_mp, const hardware::System& system,
        physics::InterfacesHandler& interfaceHandler)
{
    ::leapfrog(gm, gf, phi, phi_mp, system, interfaceHandler);
}
void physics::algorithms::leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const physics::lattices::Spinorfield_eo& phi, const physics::lattices::Spinorfield_eo& phi_mp, const hardware::System& system,
        physics::InterfacesHandler& interfaceHandler)
{
    ::leapfrog(gm, gf, phi, phi_mp, system, interfaceHandler);
}

template<class SPINORFIELD> void leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const SPINORFIELD& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler)
{
    const physics::algorithms::IntegratorParametersInterface & parametersInterface = interfaceHandler.getIntegratorParametersInterface();

    logger.trace() << "\tHMC [INT]:\tstart leapfrog...";
    //it is assumed that the new gaugefield and gaugemomentum have been set to the old ones already when this function is called the first time
    if(parametersInterface.getNumTimescales() == 1) {
        leapfrog_1ts(gm, gf, phi, system, interfaceHandler);
    } else if(parametersInterface.getNumTimescales() == 2) {
        leapfrog_2ts(gm, gf, phi, system, interfaceHandler);
    } else if(parametersInterface.getNumTimescales() == 3) {
        throw Print_Error_Message("3 timescales require mass prec.");
    } else {
        throw Print_Error_Message("\tHMC [INT]:\tMore than 3 timescales is not implemented yet. Aborting...");
    }
    logger.trace() << "\tHMC [INT]:\t...finished leapfrog";
}

template<class SPINORFIELD> void leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const SPINORFIELD& phi, const SPINORFIELD& phi_mp, const hardware::System& system, physics::InterfacesHandler& interfaceHandler)
{
    const physics::algorithms::IntegratorParametersInterface & parametersInterface = interfaceHandler.getIntegratorParametersInterface();

    logger.trace() << "\tHMC [INT]:\tstart leapfrog...";
    //it is assumed that the new gaugefield and gaugemomentum have been set to the old ones already when this function is called the first time
    if(parametersInterface.getNumTimescales() == 1 || parametersInterface.getNumTimescales() == 2) {
        throw Print_Error_Message("1 or 2 timescales cannot be used with mass prec.");
    } else if(parametersInterface.getNumTimescales() == 3) {
        leapfrog_3ts(gm, gf, phi, phi_mp, system, interfaceHandler);
    } else {
        Print_Error_Message("\tHMC [INT]:\tMore than 3 timescales is not implemented yet. Aborting...");
    }
    logger.trace() << "\tHMC [INT]:\t...finished leapfrog";
}

template<class SPINORFIELD> static void leapfrog_1ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const SPINORFIELD& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler)
{
    using namespace physics::algorithms;

    const physics::algorithms::IntegratorParametersInterface & parametersInterface = interfaceHandler.getIntegratorParametersInterface();
    const physics::AdditionalParameters& additionalParameters = interfaceHandler.getAdditionalParameters<SPINORFIELD>();
    const int n0 = parametersInterface.getIntegrationSteps(0);
    const hmc_float deltaTau0 = parametersInterface.getTau() / ((hmc_float) n0);
    const hmc_float deltaTau0_half = 0.5 * deltaTau0;

    md_update_gaugemomentum(gm, deltaTau0_half, *gf, phi, system, interfaceHandler, additionalParameters);
    for (int k = 1; k < n0; k++) {
        md_update_gaugefield(gf, *gm, deltaTau0);
        md_update_gaugemomentum(gm, deltaTau0, *gf, phi, system, interfaceHandler, additionalParameters);
    }
    md_update_gaugefield(gf, *gm, deltaTau0);
    md_update_gaugemomentum(gm, deltaTau0_half, *gf, phi, system, interfaceHandler, additionalParameters);
}

template<class SPINORFIELD> static void leapfrog_2ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const SPINORFIELD& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler)
{
    using namespace physics::algorithms;

    //this uses 2 timescales (more is not implemented yet): timescale0 for the gauge-part, timescale1 for the fermion part
    //this is done after hep-lat/0209037. See also hep-lat/0506011v2 for a more advanced version
    const physics::algorithms::IntegratorParametersInterface & parametersInterface = interfaceHandler.getIntegratorParametersInterface();
    const physics::AdditionalParameters& additionalParameters = interfaceHandler.getAdditionalParameters<SPINORFIELD>();
    const int n0 = parametersInterface.getIntegrationSteps(0);
    const int n1 = parametersInterface.getIntegrationSteps(1);
    const hmc_float deltaTau1 = parametersInterface.getTau() / ((hmc_float) n1);
    const hmc_float deltaTau0 = deltaTau1 / ((hmc_float) n0);
    const hmc_float deltaTau0_half = 0.5 * deltaTau0;
    const hmc_float deltaTau1_half = 0.5 * deltaTau1;

    //this corresponds to V_s2(deltaTau/2)
    md_update_gaugemomentum_fermion(gm, deltaTau1_half, *gf, phi, system, interfaceHandler, additionalParameters);
    //now, m steps "more" are performed for the gauge-part
    //this corresponds to [V_s1(deltaTau/2/m) V_t(deltaTau/m) V_s1(deltaTau/2/m) ]^m
    for (int l = 0; l < n0; l++) {
        if(l == 0)
            md_update_gaugemomentum_gauge(gm, deltaTau0_half, *gf, system, interfaceHandler);
        md_update_gaugefield(gf, *gm, deltaTau0);
        //one has to include the case of n1=1 here
        if(l == n0 - 1 && n1 == 1)
            md_update_gaugemomentum_gauge(gm, deltaTau0_half, *gf, system, interfaceHandler);
        else
            md_update_gaugemomentum_gauge(gm, deltaTau0, *gf, system, interfaceHandler);
    }
    for (int k = 1; k < n1; k++) {
        //this corresponds to V_s2(deltaTau)
        md_update_gaugemomentum_fermion(gm, deltaTau1, *gf, phi, system, interfaceHandler, additionalParameters);
        for (int l = 0; l < n0; l++) {
            //this corresponds to [V_s1(deltaTau/2/m) V_t(deltaTau/m) V_s1(deltaTau/2/m) ]^m
            // where the first half_step has been carried out above already
            md_update_gaugefield(gf, *gm, deltaTau0);
            //md_update_gaugemomentum_gauge(deltaTau0_half);
            if(l == n0 - 1 && k == n1 - 1)
                md_update_gaugemomentum_gauge(gm, deltaTau0_half, *gf, system, interfaceHandler);
            else
                md_update_gaugemomentum_gauge(gm, deltaTau0, *gf, system, interfaceHandler);
        }
    }
    //this corresponds to the missing V_s2(deltaTau/2)
    md_update_gaugemomentum_fermion(gm, deltaTau1_half, *gf, phi, system, interfaceHandler, additionalParameters);
}

template<class SPINORFIELD> static void leapfrog_3ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const SPINORFIELD& phi, const SPINORFIELD& phi_mp, const hardware::System& system, physics::InterfacesHandler& interfaceHandler)
{
    using namespace physics::algorithms;

    // just like with 2 timescales...
    const physics::algorithms::IntegratorParametersInterface & parametersInterface = interfaceHandler.getIntegratorParametersInterface();
    const int n0 = parametersInterface.getIntegrationSteps(0);
    const int n1 = parametersInterface.getIntegrationSteps(1);
    const int n2 = parametersInterface.getIntegrationSteps(2);
    const hmc_float deltaTau2 = parametersInterface.getTau() / ((hmc_float) n2);
    const hmc_float deltaTau1 = deltaTau2 / ((hmc_float) n1);
    const hmc_float deltaTau0 = deltaTau1 / ((hmc_float) n0);
    const hmc_float deltaTau0_half = 0.5 * deltaTau0;
    const hmc_float deltaTau1_half = 0.5 * deltaTau1;
    const hmc_float deltaTau2_half = 0.5 * deltaTau2;

    //In this case one has to call the "normal" md_update_gaugemomentum_fermion with the heavier mass
    const physics::AdditionalParameters& additionalParametersMp = interfaceHandler.getAdditionalParameters<SPINORFIELD>(true);

    md_update_gaugemomentum_detratio(gm, deltaTau2_half, *gf, phi_mp, system, interfaceHandler);
    //now, n1 steps "more" are performed for the fermion-part
    for (int l = 0; l < n1; l++) {
        if(l == 0)
            md_update_gaugemomentum_fermion(gm, deltaTau1_half, *gf, phi, system, interfaceHandler, additionalParametersMp);
        //now, n0 steps "more" are performed for the gauge-part
        for (int j = 0; j < n0; j++) {
            if(l == 0 && j == 0)
            md_update_gaugemomentum_gauge(gm, deltaTau0_half, *gf, system, interfaceHandler);
            md_update_gaugefield(gf, *gm, deltaTau0);
            if(j == n0 - 1 && l == n1 - 1 && n2 == 1)
                md_update_gaugemomentum_gauge(gm, deltaTau0_half, *gf, system, interfaceHandler);
            else
                md_update_gaugemomentum_gauge(gm, deltaTau0, *gf, system, interfaceHandler);
        }
        if(l == n1 - 1 && n2 == 1)
            md_update_gaugemomentum_fermion(gm, deltaTau1_half, *gf, phi, system, interfaceHandler, additionalParametersMp);
        else
            md_update_gaugemomentum_fermion(gm, deltaTau1, *gf, phi, system, interfaceHandler, additionalParametersMp);
    }
    //perform n2 - 1 intermediate steps
    for (int k = 1; k < n2; k++) {
        md_update_gaugemomentum_detratio(gm, deltaTau2, *gf, phi_mp, system, interfaceHandler);
        for (int l = 0; l < n1; l++) {
            for (int j = 0; j < n0; j++) {
                md_update_gaugefield(gf, *gm, deltaTau0);
                if(j == n0 - 1 && l == n1 - 1 && k == n2 - 1)
                    md_update_gaugemomentum_gauge(gm, deltaTau0_half, *gf, system, interfaceHandler);
                else
                    md_update_gaugemomentum_gauge(gm, deltaTau0, *gf, system, interfaceHandler);
            }
            if(l == n1 - 1 && k == n2 - 1)
                md_update_gaugemomentum_fermion(gm, deltaTau1_half, *gf, phi, system, interfaceHandler, additionalParametersMp);
            else
                md_update_gaugemomentum_fermion(gm, deltaTau1, *gf, phi, system, interfaceHandler, additionalParametersMp);
        }
    }
    md_update_gaugemomentum_detratio(gm, deltaTau2_half, *gf, phi_mp, system, interfaceHandler);
}

void physics::algorithms::twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const physics::lattices::Spinorfield& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler)
{
    ::twomn(gm, gf, phi, system, interfaceHandler);
}
void physics::algorithms::twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler)
{
    ::twomn(gm, gf, phi, system, interfaceHandler);
}
void physics::algorithms::twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const physics::lattices::Rooted_Staggeredfield_eo& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler)
{
    ::twomn(gm, gf, phi, system, interfaceHandler);
}
void physics::algorithms::twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const physics::lattices::Spinorfield& phi, const physics::lattices::Spinorfield& phi_mp, const hardware::System& system,
        physics::InterfacesHandler& interfaceHandler)
{
    ::twomn(gm, gf, phi, phi_mp, system, interfaceHandler);
}
void physics::algorithms::twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const physics::lattices::Spinorfield_eo& phi, const physics::lattices::Spinorfield_eo& phi_mp, const hardware::System& system,
        physics::InterfacesHandler& interfaceHandler)
{
    ::twomn(gm, gf, phi, phi_mp, system, interfaceHandler);
}

template<class SPINORFIELD> void twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi,
        const hardware::System& system, physics::InterfacesHandler& interfaceHandler)
{
    const physics::algorithms::IntegratorParametersInterface & parametersInterface = interfaceHandler.getIntegratorParametersInterface();

    logger.trace() << "\tHMC [INT]\tstarting 2MN...";
    //it is assumed that the new gaugefield and gaugemomentum have been set to the old ones already when this function is called the first time
    if(parametersInterface.getNumTimescales() == 1) {
        twomn_1ts(gm, gf, phi, system, interfaceHandler);
    } else if(parametersInterface.getNumTimescales() == 2) {
        twomn_2ts(gm, gf, phi, system, interfaceHandler);
    } else if(parametersInterface.getNumTimescales() == 3) {
        throw Print_Error_Message("3 timescales require mass prec.");
    } else {
        throw Print_Error_Message("\tHMC [INT]:\tMore than 3 timescales is not implemented yet. Aborting...");
    }
    logger.debug() << "\tHMC [INT]:\tfinished 2MN";
}

template<class SPINORFIELD> void twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi,
        const SPINORFIELD& phi_mp, const hardware::System& system, physics::InterfacesHandler& interfaceHandler)
{
    const physics::algorithms::IntegratorParametersInterface & parametersInterface = interfaceHandler.getIntegratorParametersInterface();

    logger.trace() << "\tHMC [INT]\tstarting 2MN...";
    //it is assumed that the new gaugefield and gaugemomentum have been set to the old ones already when this function is called the first time
    if(parametersInterface.getNumTimescales() == 1 || parametersInterface.getNumTimescales() == 2) {
        throw Print_Error_Message("1 or 2 timescales cannot be used with mass prec.");
    } else if(parametersInterface.getNumTimescales() == 3) {
        twomn_3ts(gm, gf, phi, phi_mp, system, interfaceHandler);
    } else {
        throw Print_Error_Message("\tHMC [INT]:\tMore than 3 timescales is not implemented yet. Aborting...");
    }
    logger.debug() << "\tHMC [INT]:\tfinished 2MN";
}

template<class SPINORFIELD> void twomn_1ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const SPINORFIELD& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler)
{
    using namespace physics::algorithms;

    const physics::algorithms::IntegratorParametersInterface & parametersInterface = interfaceHandler.getIntegratorParametersInterface();
    const physics::AdditionalParameters& additionalParameters = interfaceHandler.getAdditionalParameters<SPINORFIELD>();
    const int n0 = parametersInterface.getIntegrationSteps(0);
    const hmc_float deltaTau0 = parametersInterface.getTau() / ((hmc_float) n0);
    const hmc_float deltaTau0_half = 0.5 * deltaTau0;
    const hmc_float lambda_times_deltaTau0 = deltaTau0 * parametersInterface.getLambda(0);
    const hmc_float one_minus_2_lambda = 1. - 2. * parametersInterface.getLambda(0);
    const hmc_float one_minus_2_lambda_times_deltaTau0 = one_minus_2_lambda * deltaTau0;

    md_update_gaugemomentum(gm, lambda_times_deltaTau0, *gf, phi, system, interfaceHandler, additionalParameters);
    md_update_gaugefield(gf, *gm, deltaTau0_half);
    md_update_gaugemomentum(gm, one_minus_2_lambda_times_deltaTau0, *gf, phi, system, interfaceHandler, additionalParameters);
    md_update_gaugefield(gf, *gm, deltaTau0_half);

    for (int k = 1; k < n0; k++) {
        md_update_gaugemomentum(gm, 2. * lambda_times_deltaTau0, *gf, phi, system, interfaceHandler, additionalParameters);
        md_update_gaugefield(gf, *gm, deltaTau0_half);
        md_update_gaugemomentum(gm, one_minus_2_lambda_times_deltaTau0, *gf, phi, system, interfaceHandler, additionalParameters);
        md_update_gaugefield(gf, *gm, deltaTau0_half);
    }
    md_update_gaugemomentum(gm, lambda_times_deltaTau0, *gf, phi, system, interfaceHandler, additionalParameters);
}

template<class SPINORFIELD> void twomn_2ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const SPINORFIELD& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler)
{
    using namespace physics::algorithms;

    //this is done after hep-lat/0209037. See also hep-lat/0506011v2 for a more advanced version
    const physics::algorithms::IntegratorParametersInterface & parametersInterface = interfaceHandler.getIntegratorParametersInterface();
    const physics::AdditionalParameters& additionalParameters = interfaceHandler.getAdditionalParameters<SPINORFIELD>();
    const int n0 = parametersInterface.getIntegrationSteps(0);
    const int n1 = parametersInterface.getIntegrationSteps(1);

    //this uses 2 timescales (more is not implemented yet): timescale1 for the gauge-part, timescale2 for the fermion part
    const hmc_float deltaTau1 = parametersInterface.getTau() / ((hmc_float) n1);
    //NOTE: With 2MN, the stepsize for the lower integration step is deltaTau1/(2 N0)!!
    const hmc_float deltaTau0 = deltaTau1 / (2. * (hmc_float) n0);
    const hmc_float deltaTau0_half = 0.5 * deltaTau0;
    //hmc_float deltaTau1_half = 0.5 * deltaTau1;
    const hmc_float lambda0_times_deltaTau0 = deltaTau0 * parametersInterface.getLambda(0);
    const hmc_float lambda1_times_deltaTau1 = deltaTau1 * parametersInterface.getLambda(1);
    const hmc_float one_minus_2_lambda0 = 1. - 2. * parametersInterface.getLambda(0);
    const hmc_float one_minus_2_lambda1 = 1. - 2. * parametersInterface.getLambda(1);
    const hmc_float one_minus_2_lambda0_times_deltaTau0 = one_minus_2_lambda0 * deltaTau0;
    const hmc_float one_minus_2_lambda1_times_deltaTau1 = one_minus_2_lambda1 * deltaTau1;

    md_update_gaugemomentum_fermion(gm, lambda1_times_deltaTau1, *gf, phi, system, interfaceHandler, additionalParameters);
    //now, n0 steps "more" are performed for the gauge-part
    //this corresponds to [exp(lambda*eps T(V_gauge) ) exp( eps/2 V ) exp( (1 - 2lamdba) *eps T(V_gauge) ) exp( eps/2 V ) exp( lamdba*eps T(V_ga    uge) ) ]^m
    for (int l = 0; l < n0; l++) {
        if(l == 0)
            md_update_gaugemomentum_gauge(gm, lambda0_times_deltaTau0, *gf, system, interfaceHandler);
        md_update_gaugefield(gf, *gm, deltaTau0_half);
        md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system, interfaceHandler);
        md_update_gaugefield(gf, *gm, deltaTau0_half);
        md_update_gaugemomentum_gauge(gm, 2. * lambda0_times_deltaTau0, *gf, system, interfaceHandler);
    }
    //this corresponds to V_s2( ( 1 - 2lambda) *deltaTau)
    md_update_gaugemomentum_fermion(gm, one_minus_2_lambda1_times_deltaTau1, *gf, phi, system, interfaceHandler, additionalParameters);
    //now, m steps "more" are performed for the gauge-part (again)
    for (int l = 0; l < n0; l++) {
        //the first half step has been carried out above already
        md_update_gaugefield(gf, *gm, deltaTau0_half);
        md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system, interfaceHandler);
        md_update_gaugefield(gf, *gm, deltaTau0_half);
        //in case one does not perform intermediate steps after this, one must perform a half_step only!
        if(l == n0 - 1 && n1 == 1)
            md_update_gaugemomentum_gauge(gm, lambda0_times_deltaTau0, *gf, system, interfaceHandler);
        else
            md_update_gaugemomentum_gauge(gm, 2. * lambda0_times_deltaTau0, *gf, system, interfaceHandler);
    }
    //the last V_s2(lambda*deltaTau) can be pulled into the intermediate steps
    for (int k = 1; k < n1; k++) {
        //this corresponds to V_s2(deltaTau)
        md_update_gaugemomentum_fermion(gm, 2. * lambda1_times_deltaTau1, *gf, phi, system, interfaceHandler, additionalParameters);
        for (int l = 0; l < n0; l++) {
            //the first half step has been carried out above already
            md_update_gaugefield(gf, *gm, deltaTau0_half);
            md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system, interfaceHandler);
            md_update_gaugefield(gf, *gm, deltaTau0_half);
            md_update_gaugemomentum_gauge(gm, 2. * lambda0_times_deltaTau0, *gf, system, interfaceHandler);
        }
        md_update_gaugemomentum_fermion(gm, one_minus_2_lambda1_times_deltaTau1, *gf, phi, system, interfaceHandler, additionalParameters);
        for (int l = 0; l < n0; l++) {
            //the first half step has been carried out above already
            md_update_gaugefield(gf, *gm, deltaTau0_half);
            md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system, interfaceHandler);
            md_update_gaugefield(gf, *gm, deltaTau0_half);
            if(l == n0 - 1 && k == n1 - 1)
                md_update_gaugemomentum_gauge(gm, lambda0_times_deltaTau0, *gf, system, interfaceHandler);
            else
                md_update_gaugemomentum_gauge(gm, 2. * lambda0_times_deltaTau0, *gf, system, interfaceHandler);
        }
    }
    //this corresponds to the missing V_s2(lambda*deltaTau)
    md_update_gaugemomentum_fermion(gm, lambda1_times_deltaTau1, *gf, phi, system, interfaceHandler, additionalParameters);
}

template<class SPINORFIELD> void twomn_3ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const SPINORFIELD& phi, const SPINORFIELD& phi_mp, const hardware::System& system, physics::InterfacesHandler& interfaceHandler)
{
    using namespace physics::algorithms;

    //just like with 2 timescales...
    const physics::algorithms::IntegratorParametersInterface & parametersInterface = interfaceHandler.getIntegratorParametersInterface();
    const int n0 = parametersInterface.getIntegrationSteps(0);
    const int n1 = parametersInterface.getIntegrationSteps(1);
    const int n2 = parametersInterface.getIntegrationSteps(2);

    const hmc_float deltaTau2 = parametersInterface.getTau() / ((hmc_float) n2);
    //NOTE: With 2MN, the stepsize for the lower integration step is deltaTau1/(2 Ni-1)!!
    const hmc_float deltaTau1 = deltaTau2 / (2. * (hmc_float) n1);
    const hmc_float deltaTau0 = deltaTau1 / (2. * (hmc_float) n0);
    const hmc_float deltaTau0_half = 0.5 * deltaTau0;
    const hmc_float lambda0_times_deltaTau0 = deltaTau0 * parametersInterface.getLambda(0);
    const hmc_float lambda1_times_deltaTau1 = deltaTau1 * parametersInterface.getLambda(1);
    const hmc_float lambda2_times_deltaTau2 = deltaTau2 * parametersInterface.getLambda(2);
    const hmc_float one_minus_2_lambda0 = 1. - 2. * parametersInterface.getLambda(0);
    const hmc_float one_minus_2_lambda1 = 1. - 2. * parametersInterface.getLambda(1);
    const hmc_float one_minus_2_lambda2 = 1. - 2. * parametersInterface.getLambda(2);
    const hmc_float one_minus_2_lambda0_times_deltaTau0 = one_minus_2_lambda0 * deltaTau0;
    const hmc_float one_minus_2_lambda1_times_deltaTau1 = one_minus_2_lambda1 * deltaTau1;
    const hmc_float one_minus_2_lambda2_times_deltaTau2 = one_minus_2_lambda2 * deltaTau2;

    //In this case one has to call the "normal" md_update_gaugemomentum_fermion with the heavier mass
    const physics::AdditionalParameters& additionalParametersMp = interfaceHandler.getAdditionalParameters<SPINORFIELD>(true);

    md_update_gaugemomentum_detratio(gm, lambda2_times_deltaTau2, *gf, phi_mp, system, interfaceHandler);
    for (int l = 0; l < n1; l++) {
        if(l == 0)
            md_update_gaugemomentum_fermion(gm, lambda1_times_deltaTau1, *gf, phi, system, interfaceHandler, additionalParametersMp);
        for (int j = 0; j < n0; j++) {
            if(j == 0 && l == 0)
                md_update_gaugemomentum_gauge(gm, lambda0_times_deltaTau0, *gf, system, interfaceHandler);
            md_update_gaugefield(gf, *gm, deltaTau0_half);
            md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system, interfaceHandler);
            md_update_gaugefield(gf, *gm, deltaTau0_half);
            md_update_gaugemomentum_gauge(gm, 2. * lambda0_times_deltaTau0, *gf, system, interfaceHandler);
        }
        md_update_gaugemomentum_fermion(gm, one_minus_2_lambda1_times_deltaTau1, *gf, phi, system, interfaceHandler, additionalParametersMp);
        for (int j = 0; j < n0; j++) {
            md_update_gaugefield(gf, *gm, deltaTau0_half);
            md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system, interfaceHandler);
            md_update_gaugefield(gf, *gm, deltaTau0_half);
            md_update_gaugemomentum_gauge(gm, 2. * lambda0_times_deltaTau0, *gf, system, interfaceHandler);
        }
        md_update_gaugemomentum_fermion(gm, 2. * lambda1_times_deltaTau1, *gf, phi, system, interfaceHandler, additionalParametersMp);
    }
    md_update_gaugemomentum_detratio(gm, one_minus_2_lambda2_times_deltaTau2, *gf, phi_mp, system, interfaceHandler);
    for (int l = 0; l < n1; l++) {
        for (int j = 0; j < n0; j++) {
            md_update_gaugefield(gf, *gm, deltaTau0_half);
            md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system, interfaceHandler);
            md_update_gaugefield(gf, *gm, deltaTau0_half);
            md_update_gaugemomentum_gauge(gm, 2. * lambda0_times_deltaTau0, *gf, system, interfaceHandler);
        }
        md_update_gaugemomentum_fermion(gm, one_minus_2_lambda1_times_deltaTau1, *gf, phi, system, interfaceHandler, additionalParametersMp);
        for (int j = 0; j < n0; j++) {
            md_update_gaugefield(gf, *gm, deltaTau0_half);
            md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system, interfaceHandler);
            md_update_gaugefield(gf, *gm, deltaTau0_half);
            if(j == n0 - 1 && l == n1 - 1 && n2 == 1)
                md_update_gaugemomentum_gauge(gm, lambda0_times_deltaTau0, *gf, system, interfaceHandler);
            else
                md_update_gaugemomentum_gauge(gm, 2. * lambda0_times_deltaTau0, *gf, system, interfaceHandler);
        }
        if(l == n1 - 1 && n2 == 1)
            md_update_gaugemomentum_fermion(gm, lambda1_times_deltaTau1, *gf, phi, system, interfaceHandler, additionalParametersMp);
        else
            md_update_gaugemomentum_fermion(gm, 2. * lambda1_times_deltaTau1, *gf, phi, system, interfaceHandler, additionalParametersMp);
    }
    for (int k = 1; k < n2; k++) {
        //this corresponds to V_s2(deltaTau)
        md_update_gaugemomentum_detratio(gm, 2. * lambda2_times_deltaTau2, *gf, phi_mp, system, interfaceHandler);
        for (int l = 0; l < n1; l++) {
            for (int j = 0; j < n0; j++) {
                md_update_gaugefield(gf, *gm, deltaTau0_half);
                md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system, interfaceHandler);
                md_update_gaugefield(gf, *gm, deltaTau0_half);
                md_update_gaugemomentum_gauge(gm, 2. * lambda0_times_deltaTau0, *gf, system, interfaceHandler);
            }
            md_update_gaugemomentum_fermion(gm, one_minus_2_lambda1_times_deltaTau1, *gf, phi, system, interfaceHandler, additionalParametersMp);
            for (int j = 0; j < n0; j++) {
                md_update_gaugefield(gf, *gm, deltaTau0_half);
                md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system, interfaceHandler);
                md_update_gaugefield(gf, *gm, deltaTau0_half);
                md_update_gaugemomentum_gauge(gm, 2. * lambda0_times_deltaTau0, *gf, system, interfaceHandler);
            }
            md_update_gaugemomentum_fermion(gm, 2. * lambda1_times_deltaTau1, *gf, phi, system, interfaceHandler, additionalParametersMp);
        }
        md_update_gaugemomentum_detratio(gm, one_minus_2_lambda2_times_deltaTau2, *gf, phi_mp, system, interfaceHandler);
        for (int l = 0; l < n1; l++) {
            for (int j = 0; j < n0; j++) {
                md_update_gaugefield(gf, *gm, deltaTau0_half);
                md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system, interfaceHandler);
                md_update_gaugefield(gf, *gm, deltaTau0_half);
                md_update_gaugemomentum_gauge(gm, 2. * lambda0_times_deltaTau0, *gf, system, interfaceHandler);
            }
            md_update_gaugemomentum_fermion(gm, one_minus_2_lambda1_times_deltaTau1, *gf, phi, system, interfaceHandler, additionalParametersMp);
            for (int j = 0; j < n0; j++) {
                md_update_gaugefield(gf, *gm, deltaTau0_half);
                md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system, interfaceHandler);
                md_update_gaugefield(gf, *gm, deltaTau0_half);
                if(j == n0 - 1 && l == n1 - 1 && k == n2 - 1)
                    md_update_gaugemomentum_gauge(gm, lambda0_times_deltaTau0, *gf, system, interfaceHandler);
                else
                    md_update_gaugemomentum_gauge(gm, 2. * lambda0_times_deltaTau0, *gf, system, interfaceHandler);
            }
            if(l == n1 - 1 && k == n2 - 1)
                md_update_gaugemomentum_fermion(gm, lambda1_times_deltaTau1, *gf, phi, system, interfaceHandler, additionalParametersMp);
            else
                md_update_gaugemomentum_fermion(gm, 2. * lambda1_times_deltaTau1, *gf, phi, system, interfaceHandler, additionalParametersMp);
        }
    }
    md_update_gaugemomentum_detratio(gm, lambda2_times_deltaTau2, *gf, phi_mp, system, interfaceHandler);
}

void physics::algorithms::integrator(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const physics::lattices::Spinorfield& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler)
{
    ::integrator(gm, gf, phi, system, interfaceHandler);
}
void physics::algorithms::integrator(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const physics::lattices::Spinorfield_eo& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler)
{
    ::integrator(gm, gf, phi, system, interfaceHandler);
}
void physics::algorithms::integrator(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const physics::lattices::Rooted_Staggeredfield_eo& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler)
{
    ::integrator(gm, gf, phi, system, interfaceHandler);
}

template<class SPINORFIELD> void integrator(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const SPINORFIELD& phi, const hardware::System& system, physics::InterfacesHandler& interfaceHandler)
{
    const physics::algorithms::IntegratorParametersInterface & parametersInterface = interfaceHandler.getIntegratorParametersInterface();

    check_integrator_params(interfaceHandler);

    //CP: actual integrator calling
    switch (parametersInterface.getIntegrator(0)) {
        case common::leapfrog:
            leapfrog(gm, gf, phi, system, interfaceHandler);
            break;
        case common::twomn:
            twomn(gm, gf, phi, system, interfaceHandler);
            break;
    }
    return;
}

void physics::algorithms::integrator(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const physics::lattices::Spinorfield& phi, const physics::lattices::Spinorfield& phi_mp, const hardware::System& system,
        physics::InterfacesHandler& interfaceHandler)
{
    ::integrator(gm, gf, phi, phi_mp, system, interfaceHandler);
}
void physics::algorithms::integrator(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const physics::lattices::Spinorfield_eo& phi, const physics::lattices::Spinorfield_eo& phi_mp, const hardware::System& system,
        physics::InterfacesHandler& interfaceHandler)
{
    ::integrator(gm, gf, phi, phi_mp, system, interfaceHandler);
}

template<class SPINORFIELD> void integrator(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf,
        const SPINORFIELD& phi, const SPINORFIELD& phi_mp, const hardware::System& system, physics::InterfacesHandler& interfaceHandler)
{
    const physics::algorithms::IntegratorParametersInterface & parametersInterface = interfaceHandler.getIntegratorParametersInterface();

    check_integrator_params(interfaceHandler);

    //CP: actual integrator calling
    switch (parametersInterface.getIntegrator(0)) {
        case common::leapfrog:
            leapfrog(gm, gf, phi, phi_mp, system, interfaceHandler);
            break;
        case common::twomn:
            twomn(gm, gf, phi, phi_mp, system, interfaceHandler);
            break;
    }
    return;
}

static void check_integrator_params(physics::InterfacesHandler& interfaceHandler)
{

    const physics::algorithms::IntegratorParametersInterface & parametersInterface = interfaceHandler.getIntegratorParametersInterface();
    auto const timescales = parametersInterface.getNumTimescales();

    //CP: at the moment, one can only use the same type of integrator if one uses more then one timescale...
    if(timescales == 2) {
        if((parametersInterface.getIntegrator(0) != parametersInterface.getIntegrator(1))) {
            throw Print_Error_Message("\tHMC [INT]:\tDifferent timescales must use the same integrator up to now!\nAborting...");
        }
    }
    if(timescales == 3) {
        if((parametersInterface.getIntegrator(0) != parametersInterface.getIntegrator(1) || parametersInterface.getIntegrator(0) != parametersInterface.getIntegrator(2))) {
            throw Print_Error_Message("\tHMC [INT]:\tDifferent timescales must use the same integrator up to now!\nAborting...");
        }
    }
    //CP: check if one of the integrationsteps is 0. This would lead to a divison by zero!
    logger.debug() << "timescales = " << timescales;
    switch (timescales) {
        case 1:
            if(parametersInterface.getIntegrationSteps(0) == 0) {
                throw Print_Error_Message("\tHMC [INT]:\tNumber of integrationsteps cannot be zero! Check settings!\nAborting...");
            }
            break;
        case 2:
            if(parametersInterface.getIntegrationSteps(0) == 0 || parametersInterface.getIntegrationSteps(1) == 0) {
                throw Print_Error_Message("\tHMC [INT]:\tNumber of integrationsteps cannot be zero! Check settings!\nAborting...");
            }
            break;
        case 3:
            if(parametersInterface.getIntegrationSteps(0) == 0 || parametersInterface.getIntegrationSteps(1) == 0 || parametersInterface.getIntegrationSteps(2) == 0) {
                throw Print_Error_Message("\tHMC [INT]:\tNumber of integrationsteps cannot be zero! Check settings!\nAborting...");
            }
            break;
    }
    //CP: check if 2 ts are used with mass-preconditioning or 3 ts without mass-preconditioning. In these cases the program does not behave well defined, since this is all
    //    hardcoded
    ///@todo This will not be needed if the integration is restructured!
    if((timescales == 3 && parametersInterface.getUseMp() == false) || (timescales == 2 && parametersInterface.getUseMp() == true)) {
        throw Print_Error_Message(
                "\tHMC [INT]:\tSetting for mass-preconditioning and number of timescales do not fit!\nUse either mass-preconditioning and 3 timescales or no mass-preonditioning and 2 timescales!\nAborting...");
    }
}
