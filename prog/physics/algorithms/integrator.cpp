/** @file
 * Implementation of the integrator algorithms
 *
 * (c) 2013 Matthias Bach <bach@compeng.uni-frankfurt.de>
 * (c) 2012-2013 Christopher Pinke <pinke@compeng.uni-frankfurt.de>
 */

#include "integrator.hpp"
#include "molecular_dynamics.hpp"
#include "../../meta/util.hpp"

template<class SPINORFIELD> static void leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system);
template<class SPINORFIELD> static void leapfrog_1ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system);
template<class SPINORFIELD> static void leapfrog_2ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system);
template<class SPINORFIELD> static void leapfrog_3ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system);

void physics::algorithms::leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const physics::lattices::Spinorfield& phi, const hardware::System& system)
{
	::leapfrog(gm, gf, phi, system);
}
void physics::algorithms::leapfrog(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const physics::lattices::Spinorfield_eo& phi, const hardware::System& system)
{
	::leapfrog(gm, gf, phi, system);
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
		leapfrog_3ts(gm, gf, phi, system);
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

template<class SPINORFIELD> static void leapfrog_3ts(const physics::lattices::Gaugemomenta * const gm, const physics::lattices::Gaugefield * const gf, const SPINORFIELD& phi, const hardware::System& system)
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

	md_update_gaugemomentum_detratio(gm, deltaTau2_half, *gf, phi, system);
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
		md_update_gaugemomentum_detratio(gm, deltaTau2, *gf, phi, system);
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
	md_update_gaugemomentum_detratio(gm, deltaTau2_half, *gf, phi, system);
}
