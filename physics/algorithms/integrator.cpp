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

template<class SPINORFIELD> static void integrator(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system);
template<class SPINORFIELD> static void integrator(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const SPINORFIELD& phi_mp, const hardware::System& system);

template<class SPINORFIELD> static void leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system);
template<class SPINORFIELD> static void leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const SPINORFIELD& phi_mp, const hardware::System& system);
template<class SPINORFIELD> static void leapfrog_1ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system);
template<class SPINORFIELD> static void leapfrog_2ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system);
template<class SPINORFIELD> static void leapfrog_3ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const SPINORFIELD& phi_mp, const hardware::System& system);

template<class SPINORFIELD> static void twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system);
template<class SPINORFIELD> static void twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const SPINORFIELD& phi_mp, const hardware::System& system);
template<class SPINORFIELD> static void twomn_1ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system);
template<class SPINORFIELD> static void twomn_2ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system);
template<class SPINORFIELD> static void twomn_3ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const SPINORFIELD& phi_mp, const hardware::System& system);

static void check_integrator_params(const meta::Inputparameters& params);


void physics::algorithms::leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const physics::lattices::Spinorfield& phi, const hardware::System& system)
{
	::leapfrog(gm, gf, phi, system);
}
void physics::algorithms::leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system)
{
	::leapfrog(gm, gf, phi, system);
}
void physics::algorithms::leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const physics::lattices::Rooted_Staggeredfield_eo& phi, const hardware::System& system)
{
	::leapfrog(gm, gf, phi, system);
}
void physics::algorithms::leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const physics::lattices::Spinorfield& phi, const physics::lattices::Spinorfield& phi_mp, const hardware::System& system)
{
	::leapfrog(gm, gf, phi, phi_mp, system);
}
void physics::algorithms::leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const physics::lattices::Spinorfield_eo& phi, const physics::lattices::Spinorfield_eo& phi_mp, const hardware::System& system)
{
	::leapfrog(gm, gf, phi, phi_mp, system);
}

template<class SPINORFIELD> void leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system)
{
	auto params = system.get_inputparameters();

	logger.trace() << "\tHMC [INT]:\tstart leapfrog...";
	//it is assumed that the new gaugefield and gaugemomentum have been set to the old ones already when this function is called the first time
	if(params.get_num_timescales() == 1) {
		leapfrog_1ts(gm, gf, phi, system);
	} else if (params.get_num_timescales() == 2) {
		leapfrog_2ts(gm, gf, phi, system);
	} else if (params.get_num_timescales() == 3) {
		throw Print_Error_Message("3 timescales require mass prec.");
	} else {
		throw Print_Error_Message("\tHMC [INT]:\tMore than 3 timescales is not implemented yet. Aborting...");
	}
	logger.trace() << "\tHMC [INT]:\t...finished leapfrog";
}

template<class SPINORFIELD> void leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const SPINORFIELD& phi_mp, const hardware::System& system)
{
	auto params = system.get_inputparameters();

	logger.trace() << "\tHMC [INT]:\tstart leapfrog...";
	//it is assumed that the new gaugefield and gaugemomentum have been set to the old ones already when this function is called the first time
	if(params.get_num_timescales() == 1 || params.get_num_timescales() == 2) {
		throw Print_Error_Message("1 or 2 timescales cannot be used with mass prec.");
	} else if (params.get_num_timescales() == 3) {
		leapfrog_3ts(gm, gf, phi, phi_mp, system);
	} else {
		Print_Error_Message("\tHMC [INT]:\tMore than 3 timescales is not implemented yet. Aborting...");
	}
	logger.trace() << "\tHMC [INT]:\t...finished leapfrog";
}

template<class SPINORFIELD> static void leapfrog_1ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system)
{
	using namespace physics::algorithms;

	auto params = system.get_inputparameters();
	const int n0 = params.get_integrationsteps(0);
	const hmc_float deltaTau0 = params.get_tau() / ((hmc_float) n0);
	const hmc_float deltaTau0_half = 0.5 * deltaTau0;

	md_update_gaugemomentum(gm, deltaTau0_half, *gf, phi, system);
	for(int k = 1; k < n0; k++) {
		md_update_gaugefield(gf, *gm, deltaTau0);
		md_update_gaugemomentum(gm, deltaTau0, *gf, phi, system);
	}
	md_update_gaugefield(gf, *gm, deltaTau0);
	md_update_gaugemomentum(gm, deltaTau0_half, *gf, phi, system);
}

template<class SPINORFIELD> static void leapfrog_2ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system)
{
	using namespace physics::algorithms;

	//this uses 2 timescales (more is not implemented yet): timescale0 for the gauge-part, timescale1 for the fermion part
	//this is done after hep-lat/0209037. See also hep-lat/0506011v2 for a more advanced version
	auto params = system.get_inputparameters();
	const int n0 = params.get_integrationsteps(0);
	const int n1 = params.get_integrationsteps(1);
	const hmc_float deltaTau1 = params.get_tau() / ((hmc_float) n1);
	const hmc_float deltaTau0 = deltaTau1 / ( (hmc_float) n0 );
	const hmc_float deltaTau0_half = 0.5 * deltaTau0;
	const hmc_float deltaTau1_half = 0.5 * deltaTau1;

	//this corresponds to V_s2(deltaTau/2)
	md_update_gaugemomentum_fermion(gm, deltaTau1_half, *gf, phi, system);
	//now, m steps "more" are performed for the gauge-part
	//this corresponds to [V_s1(deltaTau/2/m) V_t(deltaTau/m) V_s1(deltaTau/2/m) ]^m
	for(int l = 0; l < n0; l++) {
		if(l == 0) md_update_gaugemomentum_gauge(gm, deltaTau0_half, *gf, system);
		md_update_gaugefield(gf, *gm, deltaTau0);
		//one has to include the case of n1=1 here
		if(l == n0 - 1 && n1 == 1) md_update_gaugemomentum_gauge(gm, deltaTau0_half, *gf, system);
		else md_update_gaugemomentum_gauge(gm, deltaTau0, *gf, system);
	}
	for(int k = 1; k < n1; k++) {
		//this corresponds to V_s2(deltaTau)
		md_update_gaugemomentum_fermion(gm, deltaTau1, *gf, phi, system);
		for(int l = 0; l < n0; l++) {
			//this corresponds to [V_s1(deltaTau/2/m) V_t(deltaTau/m) V_s1(deltaTau/2/m) ]^m
			// where the first half_step has been carried out above already
			md_update_gaugefield(gf, *gm, deltaTau0);
			//md_update_gaugemomentum_gauge(deltaTau0_half);
			if(l == n0 - 1 && k == n1 - 1) md_update_gaugemomentum_gauge(gm, deltaTau0_half, *gf, system);
			else md_update_gaugemomentum_gauge(gm, deltaTau0, *gf, system);
		}
	}
	//this corresponds to the missing V_s2(deltaTau/2)
	md_update_gaugemomentum_fermion(gm, deltaTau1_half, *gf, phi, system);
}

template<class SPINORFIELD> static void leapfrog_3ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const SPINORFIELD& phi_mp, const hardware::System& system)
{
	using namespace physics::algorithms;

	// just like with 2 timescales...
	auto params = system.get_inputparameters();
	const int n0 = params.get_integrationsteps(0);
	const int n1 = params.get_integrationsteps(1);
	const int n2 = params.get_integrationsteps(2);
	const hmc_float deltaTau2 = params.get_tau() / ((hmc_float) n2);
	const hmc_float deltaTau1 = deltaTau2 / ( (hmc_float) n1 );
	const hmc_float deltaTau0 = deltaTau1 / ( (hmc_float) n0 );
	const hmc_float deltaTau0_half = 0.5 * deltaTau0;
	const hmc_float deltaTau1_half = 0.5 * deltaTau1;
	const hmc_float deltaTau2_half = 0.5 * deltaTau2;

	//In this case one has to call the "normal" md_update_gaugemomentum_fermion with the heavier mass
	const hmc_float kappa_tmp = params.get_kappa_mp();
	const hmc_float mubar_tmp = meta::get_mubar_mp(params);

	md_update_gaugemomentum_detratio(gm, deltaTau2_half, *gf, phi_mp, system);
	//now, n1 steps "more" are performed for the fermion-part
	for(int l = 0; l < n1; l++) {
		if(l == 0) md_update_gaugemomentum_fermion(gm, deltaTau1_half, *gf, phi, system, kappa_tmp, mubar_tmp);
		//now, n0 steps "more" are performed for the gauge-part
		for(int j = 0; j < n0; j++) {
			if(l == 0 && j == 0) md_update_gaugemomentum_gauge(gm, deltaTau0_half, *gf, system);
			md_update_gaugefield(gf, *gm, deltaTau0);
			if(j == n0 - 1 && l == n1 - 1 && n2 == 1) md_update_gaugemomentum_gauge(gm, deltaTau0_half, *gf, system);
			else md_update_gaugemomentum_gauge(gm, deltaTau0, *gf, system);
		}
		if(l == n1 - 1 && n2 == 1) md_update_gaugemomentum_fermion(gm, deltaTau1_half, *gf, phi, system, kappa_tmp, mubar_tmp);
		else md_update_gaugemomentum_fermion(gm, deltaTau1, *gf, phi, system, kappa_tmp, mubar_tmp);
	}
	//perform n2 - 1 intermediate steps
	for(int k = 1; k < n2; k++) {
		md_update_gaugemomentum_detratio(gm, deltaTau2, *gf, phi_mp, system);
		for(int l = 0; l < n1; l++) {
			for(int j = 0; j < n0; j++) {
				md_update_gaugefield(gf, *gm, deltaTau0);
				if(j == n0 - 1 && l == n1 - 1 &&  k == n2 - 1) md_update_gaugemomentum_gauge(gm, deltaTau0_half, *gf, system);
				else md_update_gaugemomentum_gauge(gm, deltaTau0, *gf, system);
			}
			if(l == n1 - 1 && k == n2 - 1) md_update_gaugemomentum_fermion(gm, deltaTau1_half, *gf, phi, system, kappa_tmp, mubar_tmp);
			else md_update_gaugemomentum_fermion(gm, deltaTau1, *gf, phi, system, kappa_tmp, mubar_tmp);
		}
	}
	md_update_gaugemomentum_detratio(gm, deltaTau2_half, *gf, phi_mp, system);
}

void physics::algorithms::twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const physics::lattices::Spinorfield& phi, const hardware::System& system)
{
	::twomn(gm, gf, phi, system);
}
void physics::algorithms::twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system)
{
	::twomn(gm, gf, phi, system);
}
void physics::algorithms::twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const physics::lattices::Rooted_Staggeredfield_eo& phi, const hardware::System& system)
{
	::twomn(gm, gf, phi, system);
}
void physics::algorithms::twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const physics::lattices::Spinorfield& phi, const physics::lattices::Spinorfield& phi_mp, const hardware::System& system)
{
	::twomn(gm, gf, phi, phi_mp, system);
}
void physics::algorithms::twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const physics::lattices::Spinorfield_eo& phi, const physics::lattices::Spinorfield_eo& phi_mp, const hardware::System& system)
{
	::twomn(gm, gf, phi, phi_mp, system);
}

template<class SPINORFIELD> void twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system)
{
	auto params = system.get_inputparameters();

	logger.trace() << "\tHMC [INT]\tstarting 2MN...";
	//it is assumed that the new gaugefield and gaugemomentum have been set to the old ones already when this function is called the first time
	if(params.get_num_timescales() == 1) {
		twomn_1ts(gm, gf, phi, system);
	} else if (params.get_num_timescales() == 2) {
		twomn_2ts(gm, gf, phi, system);
	} else if (params.get_num_timescales() == 3) {
		throw Print_Error_Message("3 timescales require mass prec.");
	} else {
		throw Print_Error_Message("\tHMC [INT]:\tMore than 3 timescales is not implemented yet. Aborting...");
	}
	logger.debug() << "\tHMC [INT]:\tfinished 2MN";
}

template<class SPINORFIELD> void twomn(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const SPINORFIELD& phi_mp, const hardware::System& system)
{
	auto params = system.get_inputparameters();

	logger.trace() << "\tHMC [INT]\tstarting 2MN...";
	//it is assumed that the new gaugefield and gaugemomentum have been set to the old ones already when this function is called the first time
	if(params.get_num_timescales() == 1 || params.get_num_timescales() == 2) {
		throw Print_Error_Message("1 or 2 timescales cannot be used with mass prec.");
	} else if (params.get_num_timescales() == 3) {
		twomn_3ts(gm, gf, phi, phi_mp, system);
	} else {
		throw Print_Error_Message("\tHMC [INT]:\tMore than 3 timescales is not implemented yet. Aborting...");
	}
	logger.debug() << "\tHMC [INT]:\tfinished 2MN";
}

template<class SPINORFIELD> void twomn_1ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system)
{
	using namespace physics::algorithms;

	auto params = system.get_inputparameters();
	const int n0 = params.get_integrationsteps(0);
	const hmc_float deltaTau0 = params.get_tau() / ((hmc_float) n0);
	const hmc_float deltaTau0_half = 0.5 * deltaTau0;
	const hmc_float lambda_times_deltaTau0 = deltaTau0 * params.get_lambda(0);
	const hmc_float one_minus_2_lambda = 1. - 2.*params.get_lambda(0);
	const hmc_float one_minus_2_lambda_times_deltaTau0 = one_minus_2_lambda * deltaTau0;

	md_update_gaugemomentum(gm, lambda_times_deltaTau0, *gf, phi, system);
	md_update_gaugefield(gf, *gm, deltaTau0_half);
	md_update_gaugemomentum(gm, one_minus_2_lambda_times_deltaTau0, *gf, phi, system);
	md_update_gaugefield(gf, *gm, deltaTau0_half);

	for(int k = 1; k < n0; k++) {
		md_update_gaugemomentum(gm, 2.*lambda_times_deltaTau0, *gf, phi, system);
		md_update_gaugefield(gf, *gm, deltaTau0_half);
		md_update_gaugemomentum(gm, one_minus_2_lambda_times_deltaTau0, *gf, phi, system);
		md_update_gaugefield(gf, *gm, deltaTau0_half);
	}
	md_update_gaugemomentum(gm, lambda_times_deltaTau0, *gf, phi, system);
}

template<class SPINORFIELD> void twomn_2ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system)
{
	using namespace physics::algorithms;

	//this is done after hep-lat/0209037. See also hep-lat/0506011v2 for a more advanced version
	auto params = system.get_inputparameters();
	const int n0 = params.get_integrationsteps(0);
	const int n1 = params.get_integrationsteps(1);

	//this uses 2 timescales (more is not implemented yet): timescale1 for the gauge-part, timescale2 for the fermion part
	const hmc_float deltaTau1 = params.get_tau() / ((hmc_float) n1);
	//NOTE: With 2MN, the stepsize for the lower integration step is deltaTau1/(2 N0)!!
	const hmc_float deltaTau0 = deltaTau1 / ( 2.* (hmc_float) n0 );
	const hmc_float deltaTau0_half = 0.5 * deltaTau0;
	//hmc_float deltaTau1_half = 0.5 * deltaTau1;
	const hmc_float lambda0_times_deltaTau0 = deltaTau0 * params.get_lambda(0);
	const hmc_float lambda1_times_deltaTau1 = deltaTau1 * params.get_lambda(1);
	const hmc_float one_minus_2_lambda0 = 1. - 2.*params.get_lambda(0);
	const hmc_float one_minus_2_lambda1 = 1. - 2.*params.get_lambda(1);
	const hmc_float one_minus_2_lambda0_times_deltaTau0 = one_minus_2_lambda0 * deltaTau0;
	const hmc_float one_minus_2_lambda1_times_deltaTau1 = one_minus_2_lambda1 * deltaTau1;

	md_update_gaugemomentum_fermion(gm, lambda1_times_deltaTau1, *gf, phi, system);
	//now, n0 steps "more" are performed for the gauge-part
	//this corresponds to [exp(lambda*eps T(V_gauge) ) exp( eps/2 V ) exp( (1 - 2lamdba) *eps T(V_gauge) ) exp( eps/2 V ) exp( lamdba*eps T(V_ga    uge) ) ]^m
	for(int l = 0; l < n0; l++) {
		if(l == 0) md_update_gaugemomentum_gauge(gm, lambda0_times_deltaTau0, *gf, system);
		md_update_gaugefield(gf, *gm, deltaTau0_half);
		md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system);
		md_update_gaugefield(gf, *gm, deltaTau0_half);
		md_update_gaugemomentum_gauge(gm, 2.*lambda0_times_deltaTau0, *gf, system);
	}
	//this corresponds to V_s2( ( 1 - 2lambda) *deltaTau)
	md_update_gaugemomentum_fermion(gm, one_minus_2_lambda1_times_deltaTau1, *gf, phi, system);
	//now, m steps "more" are performed for the gauge-part (again)
	for(int l = 0; l < n0; l++) {
		//the first half step has been carried out above already
		md_update_gaugefield(gf, *gm, deltaTau0_half);
		md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system);
		md_update_gaugefield(gf, *gm, deltaTau0_half);
		//in case one does not perform intermediate steps after this, one must perform a half_step only!
		if(l == n0 - 1 && n1 == 1) md_update_gaugemomentum_gauge(gm, lambda0_times_deltaTau0, *gf, system);
		else md_update_gaugemomentum_gauge(gm, 2.*lambda0_times_deltaTau0, *gf, system);
	}
	//the last V_s2(lambda*deltaTau) can be pulled into the intermediate steps
	for(int k = 1; k < n1; k++) {
		//this corresponds to V_s2(deltaTau)
		md_update_gaugemomentum_fermion(gm, 2.*lambda1_times_deltaTau1, *gf, phi, system);
		for(int l = 0; l < n0; l++) {
			//the first half step has been carried out above already
			md_update_gaugefield(gf, *gm, deltaTau0_half);
			md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system);
			md_update_gaugefield(gf, *gm, deltaTau0_half);
			md_update_gaugemomentum_gauge(gm, 2.*lambda0_times_deltaTau0, *gf, system);
		}
		md_update_gaugemomentum_fermion(gm, one_minus_2_lambda1_times_deltaTau1, *gf, phi, system);
		for(int l = 0; l < n0; l++) {
			//the first half step has been carried out above already
			md_update_gaugefield(gf, *gm, deltaTau0_half);
			md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system);
			md_update_gaugefield(gf, *gm, deltaTau0_half);
			if(l == n0 - 1 && k == n1 - 1) md_update_gaugemomentum_gauge(gm, lambda0_times_deltaTau0, *gf, system);
			else md_update_gaugemomentum_gauge(gm, 2.*lambda0_times_deltaTau0, *gf, system);
		}
	}
	//this corresponds to the missing V_s2(lambda*deltaTau)
	md_update_gaugemomentum_fermion(gm, lambda1_times_deltaTau1, *gf, phi, system);
}

template<class SPINORFIELD> void twomn_3ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const SPINORFIELD& phi_mp, const hardware::System& system)
{
	using namespace physics::algorithms;

	//just like with 2 timescales...
	auto params = system.get_inputparameters();
	const int n0 = params.get_integrationsteps(0);
	const int n1 = params.get_integrationsteps(1);
	const int n2 = params.get_integrationsteps(2);

	const hmc_float deltaTau2 = params.get_tau() / ((hmc_float) n2);
	//NOTE: With 2MN, the stepsize for the lower integration step is deltaTau1/(2 Ni-1)!!
	const hmc_float deltaTau1 = deltaTau2 / ( 2.* (hmc_float) n1 );
	const hmc_float deltaTau0 = deltaTau1 / ( 2.* (hmc_float) n0 );
	const hmc_float deltaTau1_half = 0.5 * deltaTau1;
	const hmc_float deltaTau0_half = 0.5 * deltaTau0;
	//const hmc_float deltaTau1_half = 0.5 * deltaTau1;
	const hmc_float lambda0_times_deltaTau0 = deltaTau0 * params.get_lambda(0);
	const hmc_float lambda1_times_deltaTau1 = deltaTau1 * params.get_lambda(1);
	const hmc_float lambda2_times_deltaTau2 = deltaTau2 * params.get_lambda(2);
	const hmc_float one_minus_2_lambda0 = 1. - 2.*params.get_lambda(0);
	const hmc_float one_minus_2_lambda1 = 1. - 2.*params.get_lambda(1);
	const hmc_float one_minus_2_lambda2 = 1. - 2.*params.get_lambda(2);
	const hmc_float one_minus_2_lambda0_times_deltaTau0 = one_minus_2_lambda0 * deltaTau0;
	const hmc_float one_minus_2_lambda1_times_deltaTau1 = one_minus_2_lambda1 * deltaTau1;
	const hmc_float one_minus_2_lambda2_times_deltaTau2 = one_minus_2_lambda2 * deltaTau2;

	//In this case one has to call the "normal" md_update_gaugemomentum_fermion with the heavier mass
	const hmc_float kappa = params.get_kappa_mp();
	const hmc_float mubar = meta::get_mubar_mp(params);

	md_update_gaugemomentum_detratio(gm, lambda2_times_deltaTau2, *gf, phi_mp, system);
	for(int l = 0; l < n1; l++) {
		if(l == 0) md_update_gaugemomentum_fermion(gm, lambda1_times_deltaTau1, *gf, phi, system, kappa, mubar);
		for(int j = 0; j < n0; j++) {
			if(j == 0 && l == 0) md_update_gaugemomentum_gauge(gm, lambda0_times_deltaTau0, *gf, system);
			md_update_gaugefield(gf, *gm, deltaTau0_half);
			md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system);
			md_update_gaugefield(gf, *gm, deltaTau0_half);
			md_update_gaugemomentum_gauge(gm, 2.*lambda0_times_deltaTau0, *gf, system);
		}
		md_update_gaugemomentum_fermion(gm, one_minus_2_lambda1_times_deltaTau1, *gf, phi, system, kappa, mubar);
		for(int j = 0; j < n0; j++) {
			md_update_gaugefield(gf, *gm, deltaTau0_half);
			md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system);
			md_update_gaugefield(gf, *gm, deltaTau0_half);
			md_update_gaugemomentum_gauge(gm, 2.*lambda0_times_deltaTau0, *gf, system);
		}
		md_update_gaugemomentum_fermion(gm, 2.*lambda1_times_deltaTau1, *gf, phi, system, kappa, mubar);
	}
	md_update_gaugemomentum_detratio(gm, one_minus_2_lambda2_times_deltaTau2, *gf, phi_mp, system);
	for(int l = 0; l < n1; l++) {
		for(int j = 0; j < n0; j++) {
			md_update_gaugefield(gf, *gm, deltaTau0_half);
			md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system);
			md_update_gaugefield(gf, *gm, deltaTau0_half);
			md_update_gaugemomentum_gauge(gm, 2.*lambda0_times_deltaTau0, *gf, system);
		}
		md_update_gaugemomentum_fermion(gm, one_minus_2_lambda1_times_deltaTau1, *gf, phi, system, kappa, mubar);
		for(int j = 0; j < n0; j++) {
			md_update_gaugefield(gf, *gm, deltaTau0_half);
			md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system);
			md_update_gaugefield(gf, *gm, deltaTau0_half);
			if(j == n0 - 1 && l == n1 - 1 && n2 == 1) md_update_gaugemomentum_gauge(gm, lambda0_times_deltaTau0, *gf, system);
			else md_update_gaugemomentum_gauge(gm, 2.*lambda0_times_deltaTau0, *gf, system);
		}
		if(l == n1 - 1 && n2 == 1) md_update_gaugemomentum_fermion(gm, lambda1_times_deltaTau1, *gf, phi, system, kappa, mubar);
		else md_update_gaugemomentum_fermion(gm, 2.*lambda1_times_deltaTau1, *gf, phi, system, kappa, mubar);
	}
	for(int k = 1; k < n2; k++) {
		//this corresponds to V_s2(deltaTau)
		md_update_gaugemomentum_detratio(gm, 2.*lambda2_times_deltaTau2, *gf, phi_mp, system);
		for(int l = 0; l < n1; l++) {
			for(int j = 0; j < n0; j++) {
				md_update_gaugefield(gf, *gm, deltaTau0_half);
				md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system);
				md_update_gaugefield(gf, *gm, deltaTau0_half);
				md_update_gaugemomentum_gauge(gm, 2.*lambda0_times_deltaTau0, *gf, system);
			}
			md_update_gaugemomentum_fermion(gm, one_minus_2_lambda1_times_deltaTau1, *gf, phi, system, kappa, mubar);
			for(int j = 0; j < n0; j++) {
				md_update_gaugefield(gf, *gm, deltaTau0_half);
				md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system);
				md_update_gaugefield(gf, *gm, deltaTau0_half);
				md_update_gaugemomentum_gauge(gm, 2.*lambda0_times_deltaTau0, *gf, system);
			}
			md_update_gaugemomentum_fermion(gm, 2.*lambda1_times_deltaTau1, *gf, phi, system, kappa, mubar);
		}
		md_update_gaugemomentum_detratio(gm, one_minus_2_lambda2_times_deltaTau2, *gf, phi_mp, system);
		for(int l = 0; l < n1; l++) {
			for(int j = 0; j < n0; j++) {
				md_update_gaugefield(gf, *gm, deltaTau0_half);
				md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system);
				md_update_gaugefield(gf, *gm, deltaTau0_half);
				md_update_gaugemomentum_gauge(gm, 2.*lambda0_times_deltaTau0, *gf, system);
			}
			md_update_gaugemomentum_fermion(gm, one_minus_2_lambda1_times_deltaTau1, *gf, phi, system, kappa, mubar);
			for(int j = 0; j < n0; j++) {
				md_update_gaugefield(gf, *gm, deltaTau0_half);
				md_update_gaugemomentum_gauge(gm, one_minus_2_lambda0_times_deltaTau0, *gf, system);
				md_update_gaugefield(gf, *gm, deltaTau0_half);
				if(j == n0 - 1 && l == n1 - 1 && k == n2 - 1) md_update_gaugemomentum_gauge(gm, lambda0_times_deltaTau0, *gf, system);
				else md_update_gaugemomentum_gauge(gm, 2.*lambda0_times_deltaTau0, *gf, system);
			}
			if(l == n1 - 1 && k == n2 - 1) md_update_gaugemomentum_fermion(gm, lambda1_times_deltaTau1, *gf, phi, system, kappa, mubar);
			else md_update_gaugemomentum_fermion(gm, 2.*lambda1_times_deltaTau1, *gf, phi, system, kappa, mubar);
		}
	}
	md_update_gaugemomentum_detratio(gm, lambda2_times_deltaTau2, *gf, phi_mp, system);
}

void physics::algorithms::integrator(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const physics::lattices::Spinorfield& phi, const hardware::System& system)
{
	::integrator(gm, gf, phi, system);
}
void physics::algorithms::integrator(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system)
{
	::integrator(gm, gf, phi, system);
}
void physics::algorithms::integrator(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const physics::lattices::Rooted_Staggeredfield_eo& phi, const hardware::System& system)
{
	::integrator(gm, gf, phi, system);
}

template<class SPINORFIELD> void integrator(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system)
{
	auto const params = system.get_inputparameters();

	check_integrator_params(params);

	//CP: actual integrator calling
	switch(params.get_integrator(0)) {
		case meta::Inputparameters::leapfrog:
			leapfrog(gm, gf, phi, system);
			break;
		case meta::Inputparameters::twomn:
			twomn(gm, gf, phi, system);
			break;
	}
	return;
}

void physics::algorithms::integrator(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const physics::lattices::Spinorfield& phi, const physics::lattices::Spinorfield& phi_mp, const hardware::System& system)
{
	::integrator(gm, gf, phi, phi_mp, system);
}
void physics::algorithms::integrator(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const physics::lattices::Spinorfield_eo& phi, const physics::lattices::Spinorfield_eo& phi_mp, const hardware::System& system)
{
	::integrator(gm, gf, phi, phi_mp, system);
}

template<class SPINORFIELD> void integrator(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const SPINORFIELD& phi_mp, const hardware::System& system)
{
	auto const params = system.get_inputparameters();

	check_integrator_params(params);

	//CP: actual integrator calling
	switch(params.get_integrator(0)) {
		case meta::Inputparameters::leapfrog:
			leapfrog(gm, gf, phi, phi_mp, system);
			break;
		case meta::Inputparameters::twomn:
			twomn(gm, gf, phi, phi_mp, system);
			break;
	}
	return;
}

static void check_integrator_params(const meta::Inputparameters& params)
{
	auto const timescales = params.get_num_timescales();

	//CP: at the moment, one can only use the same type of integrator if one uses more then one timescale...
	if (timescales == 2) {
		if(( params.get_integrator(0) != params.get_integrator(1)  )) {
			throw Print_Error_Message("\tHMC [INT]:\tDifferent timescales must use the same integrator up to now!\nAborting...");
		}
	}
	if (timescales == 3) {
		if(( params.get_integrator(0) != params.get_integrator(1) || params.get_integrator(0) != params.get_integrator(2) )) {
			throw Print_Error_Message("\tHMC [INT]:\tDifferent timescales must use the same integrator up to now!\nAborting...");
		}
	}
	//CP: check if one of the integrationsteps is 0. This would lead to a divison by zero!
	logger.debug() << "timescales = " << timescales;
	switch(timescales ) {
		case 1:
			if( params.get_integrationsteps(0) == 0 ) {
				throw Print_Error_Message("\tHMC [INT]:\tNumber of integrationsteps cannot be zero! Check settings!\nAborting...");
			}
			break;
		case 2:
			if( params.get_integrationsteps(0) == 0 || params.get_integrationsteps(1) == 0) {
				throw Print_Error_Message("\tHMC [INT]:\tNumber of integrationsteps cannot be zero! Check settings!\nAborting...");
			}
			break;
		case 3:
			if( params.get_integrationsteps(0) == 0 || params.get_integrationsteps(1) == 0 || params.get_integrationsteps(2) == 0 ) {
				throw Print_Error_Message("\tHMC [INT]:\tNumber of integrationsteps cannot be zero! Check settings!\nAborting...");
			}
			break;
	}
	//CP: check if 2 ts are used with mass-preconditioning or 3 ts without mass-preconditioning. In these cases the program does not behave well defined, since this is all
	//    hardcoded
	///@todo This will not be needed if the integration is restructured!
	if ( (  timescales == 3 && params.get_use_mp() == false  ) ||
	     (  timescales == 2 && params.get_use_mp() == true  ) ) {
		throw Print_Error_Message("\tHMC [INT]:\tSetting for mass-preconditioning and number of timescales do not fit!\nUse either mass-preconditioning and 3 timescales or no mass-preonditioning and 2 timescales!\nAborting...");
	}
}
