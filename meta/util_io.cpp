/** @file
 * IO utility functions
 *
 * Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 * Copyright (c) 2012 Christopher Pinke <pinke@compeng.uni-frankfurt.de>
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

#include <stdexcept>
#include "gitcommitid.h"
#include "util.hpp"
#include "../host_functionality/logger.hpp"
#include "../executables/exceptions.h"

using namespace std;

static void print_info_gauge(std::ostream* os, const meta::Inputparameters& params);
static void print_info_gauge(const meta::Inputparameters& params);
static void print_info_integrator(int number, const meta::Inputparameters& params);
static void print_info_integrator(std::ostream* os, int number, const meta::Inputparameters& params);
static void print_info_observables_fermion_io(const meta::Inputparameters& params);
static void print_info_observables_fermion_io(std::ostream * os, const meta::Inputparameters& params);
static void print_info_observables_hmc_io(const meta::Inputparameters& params);
static void print_info_observables_hmc_io(std::ostream * os, const meta::Inputparameters& params);
static void print_info_observables_rhmc_io(const meta::Inputparameters& params);
static void print_info_observables_rhmc_io(std::ostream * os, const meta::Inputparameters& params);
static void print_info_source(const meta::Inputparameters& params);
static void print_info_source(std::ostream * os, const meta::Inputparameters& params);

void meta::print_info_global(const meta::Inputparameters& params)
{
	using namespace meta;

	logger.info() << "## Build based on commit: " << GIT_COMMIT_ID;
	logger.info() << "## **********************************************************";
	logger.info() << "## Global parameters:";
	logger.info() << "## NSPACE:  " << params.get_nspace();
	logger.info() << "## NTIME:   " << params.get_ntime();
	logger.info() << "## NDIM:    " << NDIM;
	logger.info() << "## NCOLOR:  " << NC;
	if(params.get_fermact() != common::action::rooted_stagg) logger.info() << "## NSPIN:   " << NSPIN;
	logger.info() << "## **********************************************************";
	logger.info() << "## Computational parameters:";
	logger.info() << "## PREC:    " << params.get_precision();
	if(params.get_use_rec12() == true) {
		logger.info() << "## REC12:   ON";
	} else {
		logger.info() << "## REC12:   OFF";
	}
	if(params.get_use_gpu() == true) {
		logger.info() << "## USE GPU: ON";
	} else {
		logger.info() << "## USE GPU: OFF";
	}
	if(params.get_use_aniso() == true) {
		logger.info() << "## USE ANISOTROPY: ON";
	} else {
		logger.info() << "## USE ANISOTROPY: OFF";
	}
	logger.info() << "## Number of devices demanded for calculations: " << params.get_device_count()  ;
	logger.info() << "## PRNG SEED:\t" << params.get_host_seed();
	if (! ( params.get_initial_prng_state() == "" ) )
		logger.info() << "## INITIAL PRNG STATE:\t" << params.get_initial_prng_state();
	logger.info() << "## **********************************************************";
	logger.info() << "## I/O parameters:";
	logger.info() << "## Writefrequency:  " << params.get_writefrequency();
	logger.info() << "## Savefrequency:  " << params.get_savefrequency();
	switch(params.get_startcondition()) {
		case common::start_from_source: {
			std::string sf = params.get_sourcefile();
			logger.info() << "## sourcefile = " << sf;
		}
		break;
		case common::cold_start:
			logger.info() << "## COLD start";
			break;
		case common::hot_start:
			logger.info() << "## HOT start";
			break;
	}
	if(params.get_use_smearing() == 1) {
		logger.info() << "## **********************************************************";
		logger.info() << "## Apply Smearing with:";
		logger.info() << "## rho:      " << params.get_rho();
		logger.info() << "## rho_iter: " << params.get_rho_iter();
	}
}

void meta::print_info_global(std::ostream* os, const meta::Inputparameters& params)
{
	using namespace meta;

	*os  << "## Build based on commit: " << GIT_COMMIT_ID << endl;
	*os  << "## **********************************************************" << endl;
	*os  << "## Global parameters:" << endl;
	*os  << "## NSPACE:  " << params.get_nspace() << endl;
	*os  << "## NTIME:   " << params.get_ntime() << endl;
	*os  << "## NDIM:    " << NDIM << endl;
	*os  << "## NCOLOR:  " << NC << endl;
	*os  << "## NSPIN:   " << NSPIN << endl;
	*os  << "## **********************************************************" << endl;
	*os  << "## Computational parameters:" << endl;
	*os  << "## PREC:    " << params.get_precision() << endl;
	if(params.get_use_rec12() == true) {
		*os << "## REC12:   ON"  << endl;
	} else {
		*os << "## REC12:   OFF"  << endl;
	}
	if(params.get_use_gpu() == true) {
		*os << "## USE GPU: ON"  << endl;
	} else {
		*os << "## USE GPU: OFF"  << endl;
	}
	if(params.get_use_aniso() == true) {
		*os << "## USE ANISOTROPY: ON" << endl;
	} else {
		*os << "## USE ANISOTROPY: OFF" << endl;
	}
	*os  << "## Number of devices demanded for calculations: " << params.get_device_count()  << endl;
	*os  << "## PRNG SEED:\t" << params.get_host_seed() << endl;
	if (! ( params.get_initial_prng_state() == "" ) )
		*os << "## INITIAL PRNG STATE:\t" << params.get_initial_prng_state() << endl;
	*os  << "## **********************************************************" << endl;
	*os  << "## I/O parameters:" << endl;
	*os  << "## Writefrequency:  " << params.get_writefrequency() << endl;
	*os  << "## Savefrequency:  " << params.get_savefrequency() << endl;
	switch(params.get_startcondition()) {
		case common::start_from_source: {
			std::string sf = params.get_sourcefile();
			*os << "## sourcefile = " << sf << endl;;
		}
		break;
		case common::cold_start:
			*os << "## cold start" << endl;;
			break;
		case common::hot_start:
			*os << "## hot start" << endl;;
			break;
	}
	if(params.get_use_smearing() == true) {
		*os  << "## **********************************************************" << endl;
		*os  << "## Apply Smearing with:" << endl;
		*os  << "## rho:      " << params.get_rho() << endl;
		*os  << "## rho_iter: " << params.get_rho_iter() << endl;
	}
}


void meta::print_info_heatbath(const Inputparameters& params)
{
	logger.info() << "## **********************************************************";
	logger.info() << "## Simulation parameters:";
	logger.info() << "## beta           = " << params.get_beta();
	logger.info() << "## xi             = " << params.get_xi();
	logger.info() << "## thermsteps     = " << params.get_thermalizationsteps() ;
	logger.info() << "## heatbathsteps  = " << params.get_heatbathsteps();
	logger.info() << "## overrelaxsteps = " << params.get_overrelaxsteps();
	logger.info() << "## **********************************************************";
	return;
}

void meta::print_info_heatbath(std::ostream* os, const Inputparameters& params)
{
	*os  << "## **********************************************************" << endl;
	*os  << "## Simulation parameters:" << endl;
	*os  << "## beta           = " << params.get_beta() << endl;
	*os  << "## xi             = " << params.get_xi() << endl;
	*os  << "## thermsteps     = " << params.get_thermalizationsteps() << endl;
	*os  << "## heatbathsteps  = " << params.get_heatbathsteps() << endl;
	*os  << "## overrelaxsteps = " << params.get_overrelaxsteps() << endl;
	*os  << "## **********************************************************" << endl;
	return;
}

static void print_info_fermion(const meta::Inputparameters& params)
{
	using namespace meta;
	logger.info() << "## **********************************************************";
	logger.info() << "## Fermionic parameters:";
	logger.info() << "##" ;
	logger.info() << "## Boundary Conditions:";
	logger.info() << "## theta_fermion_spatial  = " << params.get_theta_fermion_spatial();
	logger.info() << "## theta_fermion_temporal = " << params.get_theta_fermion_temporal();
	logger.info() << "##" ;
	logger.info() << "## Chemical Potential:" ;
	if(params.get_use_chem_pot_re() == true)
		logger.info() << "## chem_pot_re  = " << params.get_chem_pot_re();
	else
		logger.info() << "## do not use real chem. pot.";
	if(params.get_use_chem_pot_im() == true)
		logger.info() << "## chem_pot_im = " << params.get_chem_pot_im();
	else
		logger.info() << "## do not use imag. chem. pot.";
	logger.info() << "##" ;
	if(params.get_fermact() == common::action::wilson) {
		logger.info() <<  "## fermion action: unimproved Wilson";
		logger.info() << "## kappa  = " << params.get_kappa();
	}
	if(params.get_fermact() == common::action::twistedmass) {
		logger.info() <<  "## fermion action: twisted mass Wilson";
		logger.info() << "## kappa  = " << params.get_kappa();
		logger.info() << "## mu     = " << params.get_mu();
	}
	if(params.get_fermact() == common::action::clover) {
		logger.info() <<  "## fermion action: clover Wilson";
		logger.info() << "## kappa  = " << params.get_kappa();
		logger.info() << "## csw    = " << params.get_csw();
	}
	if(params.get_fermact() == common::action::rooted_stagg) {
		logger.info() <<  "## fermion action: staggered standard ";
		logger.info() << "## mass           = " << params.get_mass();
		logger.info() << "## Num. of tastes = " << params.get_num_tastes();
	}
	logger.info() << "##" ;
	logger.info() << "## Inverter parameters:";
	logger.info() << "## precision for inversions = " << params.get_solver_prec();
	if(params.get_use_eo() == true)
		logger.info() << "## Use even-odd preconditioning" ;
	if(params.get_use_eo() == false)
		logger.info() << "## Do NOT use even-odd preconditioning";
	if(params.get_fermact() != common::action::rooted_stagg) {
		switch(params.get_solver()) {
			case common::cg:
				logger.info() << "## Use CG-solver for inversions" ;
				break;
			case common::bicgstab:
				logger.info() << "## Use BiCGStab for inversions";
				break;
			case common::bicgstab_save:
				logger.info() << "## Use BiCGStab-SAVE for inversions";
				break;
		}
		logger.info() << "## cgmax  = " << params.get_cgmax();
		logger.info() << "## iter_refresh  = " << params.get_iter_refresh();
	} else {
		logger.info() << "## Use Multi-shifted CG-solver for inversions" ;
		logger.info() << "## cgm_max  = " << params.get_cgmax();
	}

	if(params.get_profile_solver() == true)
		logger.warn() << "## Profiling of solver activated. This may influence the overall performance time!";

	if(params.get_use_merge_kernels_fermion() == true)
		logger.info() << "## Use merged fermionmatrix kernels where implemented!!";

	if(params.get_use_merge_kernels_spinor() == true)
		logger.info() << "## Use merged spinor kernels where implemented!!";

	//print extra warning if BC are set to default since this is a serious source of errors...
	if ( params.get_theta_fermion_spatial() == 0. && params.get_theta_fermion_temporal() == 0.) {
		logger.warn() << "NOTE: BCs have been set to periodic values by default!!\n\t\t\t  To change this use e.g. ThetaT/ThetaS in the input-file.";
	}
}

static void print_info_fermion(std::ostream * os, const meta::Inputparameters& params)
{
	using namespace meta;

	*os  << "## **********************************************************" << endl;
	*os  << "## Fermionic parameters:" << endl;
	*os  << "##" << endl;
	*os  << "## Boundary Conditions:" << endl;
	*os  << "## theta_fermion_spatial  = " << params.get_theta_fermion_spatial() << endl;
	*os  << "## theta_fermion_temporal = " << params.get_theta_fermion_temporal() << endl;

	*os  << "##" << endl;
	*os  << "## Chemical Potential:" << endl;
	if(params.get_use_chem_pot_re() == 1)
		*os  << "## chem_pot_re  = " << params.get_chem_pot_re() << endl;
	else
		*os  << "## do not use real chem. pot." << endl;
	if(params.get_use_chem_pot_im() == 1)
		*os  << "## chem_pot_im = " << params.get_chem_pot_im() << endl;
	else
		*os  << "## do not use imag. chem. pot." << endl;
	*os  << "##" << endl;
	if(params.get_fermact() == common::action::wilson) {
		*os <<  "## fermion action: unimproved Wilson" << endl;
		*os  << "## kappa  = " << params.get_kappa() << endl;
	}
	if(params.get_fermact() == common::action::twistedmass) {
		*os <<  "## fermion action: twisted mass Wilson" << endl;
		*os  << "## kappa  = " << params.get_kappa() << endl;
		*os  << "## mu     = " << params.get_mu() << endl;
	}
	if(params.get_fermact() == common::action::clover) {
		*os <<  "## fermion action: clover Wilson" << endl;
		*os  << "## kappa  = " << params.get_kappa() << endl;
		*os  << "## csw    = " << params.get_csw() << endl;
	}
	if(params.get_fermact() == common::action::rooted_stagg) {
		*os <<  "## fermion action: staggered standard " << endl;
		*os << "## mass           = " << params.get_mass() << endl;
		*os << "## Num. of tastes = " << params.get_num_tastes() << endl;
	}
	*os  << "##" << endl;
	*os  << "## Inverter parameters:" << endl;
	*os << "## precision for inversions = " << params.get_solver_prec() << endl;
	if(params.get_use_eo() == true)
		*os  << "## Use even-odd preconditioning" << endl;
	if(params.get_use_eo() == false)
		*os  << "## Do NOT use even-odd preconditioning" << endl;
	if(params.get_fermact() != common::action::rooted_stagg) {
		switch(params.get_solver()) {
			case common::cg:
				*os << "## Use CG-solver for inversions" << endl;
				break;
			case common::bicgstab:
				*os << "## Use BiCGStab for inversions" << endl;
				break;
			case common::bicgstab_save:
				*os << "## Use BiCGStab-SAVE for inversions" << endl;
				break;
		}
		*os << "## cgmax  = " << params.get_cgmax() << endl;
		*os << "## iter_refresh  = " << params.get_iter_refresh() << endl;
	} else {
		*os << "## Use Multi-shifted CG-solver for inversions" << endl;
		*os << "## cgm_max  = " << params.get_cgmax() << endl;
	}

	if(params.get_profile_solver() == true)
		*os << "## Profiling of solver activated. This may influence the overall performance time!" << endl;

	if(params.get_use_merge_kernels_fermion() == true)
		*os << "## Use merged fermionmatrix kernels where implemented!!" << endl;

	if(params.get_use_merge_kernels_spinor() == true)
		*os << "## Use merged spinor kernels where implemented!!" << endl;


	//print extra warning if BC are set to default since this is a serious source of errors...
	if ( params.get_theta_fermion_spatial() == 0. && params.get_theta_fermion_temporal() == 0.) {
		*os << "NOTE: BCs have been set to periodic values by default!!\n\t\t\t  To change this use e.g. ThetaT/ThetaS in the input-file." << endl;
	}
}

static void print_info_gauge(std::ostream* os, const meta::Inputparameters& params)
{
	using namespace meta;

	*os << "## **********************************************************" << endl;
	*os << "## Gauge parameters:" << endl;
	*os << "##" << endl;
	*os << "## beta:  " << params.get_beta() << endl;
	if(params.get_gaugeact() == common::action::wilson) {
		*os <<  "## gauge action: unimproved Wilson" << endl;
	}
	if(params.get_gaugeact() == common::action::tlsym) {
		*os <<  "## gauge action: tree level Symanzik" << endl;
		*os << "## c0  = " << get_c0(params) << endl;
		*os << "## c1  = " << get_c1(params) << endl;
	}
}

static void print_info_gauge(const meta::Inputparameters& params)
{
	using namespace meta;

	logger.info() << "## **********************************************************";
	logger.info() << "## Gauge parameters:";
	logger.info() << "##" ;
	logger.info() << "## beta:  " << params.get_beta();
	if(params.get_gaugeact() == common::action::wilson) {
		logger.info() <<  "## gauge action: unimproved Wilson";
	}
	if(params.get_gaugeact() == common::action::tlsym) {
		logger.info() <<  "## gauge action: tree level Symanzik";
		logger.info() << "## c0  = " << get_c0(params);
		logger.info() << "## c1  = " << get_c1(params);
	}
}

void meta::print_info_inverter(const Inputparameters& params)
{
	print_info_observables_fermion_io(params);
	print_info_fermion(params);
	print_info_source(params);
	logger.info() << "## **********************************************************";
	return;
}

void meta::print_info_inverter(std::ostream* os, const Inputparameters& params)
{
	print_info_observables_fermion_io(os, params);
	print_info_fermion(os, params);
	print_info_source(os, params);
	*os << "## **********************************************************" << endl;
	return;
}

static void print_info_integrator(int number, const meta::Inputparameters& params)
{
	using namespace meta;

	string integrator_name;
	bool print_lambda = false;
	switch(params.get_integrator(number)) {
		case common::leapfrog:
			integrator_name = "LEAPFROG";
			break;
		case common::twomn:
			integrator_name = "2MN";
			print_lambda = true;
			break;
		default:
			logger.fatal() << "Fail in getting integrator information!";
			logger.fatal() << "Aborting...";
			throw out_of_range("Unknown integrator");
	}
	logger.info() << "## integrator" << number << " = " << integrator_name;
	logger.info() << "## integrationsteps" << number << " = " << params.get_integrationsteps(number);
	if(print_lambda) logger.info() << "## lambda" << number << " = " << params.get_lambda(number);
}

static void print_info_integrator(std::ostream* os, int number, const meta::Inputparameters& params)
{
	using namespace meta;

	string integrator_name;
	bool print_lambda = false;
	switch(params.get_integrator(number)) {
		case common::leapfrog:
			integrator_name = "LEAPFROG";
			break;
		case common::twomn:
			integrator_name = "2MN";
			print_lambda = true;
			break;
		default:
			logger.fatal() << "Fail in getting integrator information!";
			logger.fatal() << "Aborting...";
			throw out_of_range("Unknown integrator");
	}
	*os << "## integrator" << number << " = " << integrator_name << endl;
	*os << "## integrationsteps" << number << " = " << params.get_integrationsteps(number) << endl;
	if(print_lambda) *os << "## lambda" << number << " = " << params.get_lambda(number) << endl;
}

void meta::print_info_hmc(const Inputparameters& params)
{
	print_info_gauge(params);
	print_info_fermion(params);
	print_info_observables_fermion_io(params);
	print_info_source(params);
	print_info_observables_hmc_io(params);
	logger.info() << "## **********************************************************";
	logger.info() << "## HMC parameters: " ;
	logger.info() << "##  ";
	logger.info() << "## tau  = " << params.get_tau();
	logger.info() << "## HMC steps  = " << params.get_hmcsteps();
	logger.info() << "## precision used in HMC-inversions = " << params.get_force_prec();
	logger.info() << "##  ";
	logger.info() << "## # Timescales  = " << params.get_num_timescales();
	if(params.get_use_gauge_only() && params.get_num_timescales() == 1) {
		logger.info() << "## use PureGaugeTheory in HMC!";
	} else if (params.get_use_gauge_only() && params.get_num_timescales() > 1) {
		logger.fatal() << "PureGaugeTheory can only be used with one timescale!\nPlease change the input-file!\nAborting...";
		throw out_of_range("Too many timescales");
	}
	//integrator infos
	for(int i = 0; i < params.get_num_timescales(); i++) {
		print_info_integrator(i, params);
	}
	if(params.get_use_mp() == true) {
		logger.info() << "##  ";
		logger.info() <<  "## use mass preconditioning:";
		if(params.get_fermact_mp() == common::action::wilson) {
			logger.info() <<  "## mp action: unimproved Wilson";
			logger.info() << "## kappa_mp  = " << params.get_kappa_mp();
		}
		if(params.get_fermact_mp() == common::action::twistedmass) {
			logger.info() <<  "## mp action: twisted mass Wilson";
			logger.info() << "## kappa_mp  = " << params.get_kappa_mp();
			logger.info() << "## mu_mp     = " << params.get_mu_mp();
		}
		if(params.get_fermact_mp() == common::action::clover) {
			logger.info() <<  "## mp action: clover Wilson";
			logger.info() << "## kappa_mp  = " << params.get_kappa_mp();
			logger.info() << "## csw_mp   = " << params.get_csw_mp();
		}
		logger.info() << "##" ;
		switch(params.get_solver_mp()) {
			case common::cg:
				logger.info() << "## Use CG-solver for mp inversions" ;
				break;
			case common::bicgstab:
				logger.info() << "## Use BiCGStab for mp inversions";
				break;
			case common::bicgstab_save:
				logger.info() << "## Use BiCGStab-SAVE for mp inversions";
				break;
		}
		logger.info() << "## cgmax_mp  = " << params.get_cgmax_mp();
		logger.info() << "## iter_refresh_mp  = " << params.get_iter_refresh_mp();
		logger.info() << "##" ;
	}
	logger.info() << "## **********************************************************";
	return;
}

void meta::print_info_hmc(std::ostream* os, const Inputparameters& params)
{
	print_info_gauge(os, params);
	print_info_fermion(os, params);
	print_info_observables_fermion_io(os, params);
	print_info_source(os, params);
	print_info_observables_hmc_io(os, params);
	*os << "## **********************************************************" << endl;
	*os << "## HMC parameters: "  << '\n';
	*os << "##  " << '\n';
	*os << "## tau  = " << params.get_tau() << '\n';
	*os << "## HMC steps  = " << params.get_hmcsteps() << '\n';
	*os << "## precision used HMC-inversions = " << params.get_force_prec() << '\n';
	*os << "##  " << '\n';
	*os << "## # Timescales  = " << params.get_num_timescales() << '\n';
	if(params.get_use_gauge_only() && params.get_num_timescales() == 1) {
		*os << "## use PureGaugeTheory in HMC!" << endl;
	} else if (params.get_use_gauge_only() && params.get_num_timescales() > 1) {
		*os << "PureGaugeTheory can only be used with one timescale!\nPlease change the input-file!\nAborting..." << endl;
		throw out_of_range("Too many timescales");
	}
	//integrator infos
	for(int i = 0; i < params.get_num_timescales(); i++) {
		print_info_integrator(os, i, params);
	}
	if(params.get_use_mp() == true) {
		*os << "##  " << '\n';
		*os <<  "## use mass preconditioning:"  << '\n';
		if(params.get_fermact_mp() == common::action::wilson) {
			*os <<  "## mp action: unimproved Wilson"  << '\n';
			*os << "## kappa_mp  = " << params.get_kappa_mp()  << '\n';
		}
		if(params.get_fermact_mp() == common::action::twistedmass) {
			*os <<  "## mp action: twisted mass Wilson"  << '\n';
			*os << "## kappa_mp  = " << params.get_kappa_mp()  << '\n';
			*os << "## mu_mp     = " << params.get_mu_mp()  << '\n';
		}
		if(params.get_fermact_mp() == common::action::clover) {
			*os <<  "## mp action: clover Wilson" << endl;
			*os << "## kappa_mp  = " << params.get_kappa_mp()  << '\n';
			*os << "## csw_mp   = " << params.get_csw_mp()  << '\n';
		}
		*os << "##"  << endl;
		switch(params.get_solver_mp()) {
			case common::cg:
				*os << "## Use CG-solver for mp inversions" << endl;
				break;
			case common::bicgstab:
				*os << "## Use BiCGStab for mp inversions" << endl;
				break;
			case common::bicgstab_save:
				*os << "## Use BiCGStab-SAVE for mp inversions" << endl;
				break;
		}
		*os << "## cgmax_mp  = " << params.get_cgmax_mp()  << '\n';
		*os << "## iter_refresh_mp  = " << params.get_iter_refresh_mp()  << '\n';
		*os << "##"  << '\n';
	}
	*os << "## **********************************************************" << '\n';
	return;
}

void meta::print_info_rhmc(const Inputparameters& params)
{
	print_info_global(params);
	print_info_configs_io(params);
	print_info_prng_io(params);
	print_info_observables_rhmc_io(params);
	print_info_gauge(params);
	print_info_fermion(params);
	logger.info() << "## **********************************************************";
	logger.info() << "## RHMC parameters: " ;
	logger.info() << "##  ";
	logger.info() << "## Rational Approximations info:";
	logger.info() << "##   - Generation of phi:";
	logger.info() << "##       + x^(+" << params.get_num_tastes() << "/8)";
	logger.info() << "##       + order = " << params.get_metro_approx_ord();
	logger.info() << "##       + range = [" << params.get_approx_lower() << " , "
	              << params.get_approx_upper() << "]";
	logger.info() << "##   - Molecular Dynamics:";
	logger.info() << "##       + x^(-" << params.get_num_tastes() << "/4)";
	logger.info() << "##       + order = " << params.get_md_approx_ord();
	logger.info() << "##       + range = [" << params.get_approx_lower() << " , "
	              << params.get_approx_upper() << "]";
	logger.info() << "##   - Evaluation of new action in Metropolis test:";
	logger.info() << "##       + x^(-" << params.get_num_tastes() << "/4)";
	logger.info() << "##       + order = " << params.get_metro_approx_ord();
	logger.info() << "##       + range = [" << params.get_approx_lower() << " , "
	              << params.get_approx_upper() << "]";
	logger.info() << "##  ";
	logger.info() << "## Strategy for finding max and min MdagM eigenvalue: " << (params.get_conservative() ? "conservative" : "NOT conservative");
	logger.info() << "##  ";
	logger.info() << "## Simulation info:";
	logger.info() << "##   - RHMC steps  = " << params.get_rhmcsteps();
	logger.info() << "##   - precision used in  Mol. Dyn.-inversions = " << params.get_force_prec();
	logger.info() << "##   - precision used in Metropolis-inversions = " << params.get_solver_prec();
	logger.info() << "##  ";
	logger.info() << "## MD integration info:";
	logger.info() << "##   - Time of integration  = " << params.get_tau();
	logger.info() << "##   - Number of Timescales  = " << params.get_num_timescales();
	if(params.get_use_gauge_only() && params.get_num_timescales() == 1) {
		logger.warn() << "## use PureGaugeTheory in RHMC!";
	} else if (params.get_use_gauge_only() && params.get_num_timescales() > 1) {
		logger.fatal() << "PureGaugeTheory can only be used with one timescale!\nPlease change the input-file!\nAborting...";
		throw out_of_range("Too many timescales");
	}
	//integrator infos
	logger.info() << "##  ";
	logger.info() << "## Integrators info:";
	for(int i = 0; i < params.get_num_timescales(); i++) {
		print_info_integrator(i, params);
	}
	if(params.get_use_mp() == true) {
		//Not yet implemented!
		/*
		logger.info() << "##  ";
		logger.info() <<  "## use mass preconditioning:";
		if(params.get_fermact_mp() == Inputparameters::wilson) {
		  logger.info() <<  "## mp action: unimproved Wilson";
		  logger.info() << "## kappa_mp  = " << params.get_kappa_mp();
		}
		if(params.get_fermact_mp() == Inputparameters::twistedmass) {
		  logger.info() <<  "## mp action: twisted mass Wilson";
		  logger.info() << "## kappa_mp  = " << params.get_kappa_mp();
		  logger.info() << "## mu_mp     = " << params.get_mu_mp();
		}
		if(params.get_fermact_mp() == Inputparameters::clover) {
		  logger.info() <<  "## mp action: clover Wilson";
		  logger.info() << "## kappa_mp  = " << params.get_kappa_mp();
		  logger.info() << "## csw_mp   = " << params.get_csw_mp();
		}
		logger.info() << "##" ;
		switch(params.get_solver_mp()) {
		  case Inputparameters::cg:
		    logger.info() << "## Use CG-solver for mp inversions" ;
		    break;
		  case Inputparameters::bicgstab:
		    logger.info() << "## Use BiCGStab for mp inversions";
		    break;
		  case Inputparameters::bicgstab_save:
		    logger.info() << "## Use BiCGStab-SAVE for mp inversions";
		    break;
		}
		logger.info() << "## cgmax_mp  = " << params.get_cgmax_mp();
		logger.info() << "## iter_refresh_mp  = " << params.get_iter_refresh_mp();
		logger.info() << "##" ;
		*/
	}
	logger.info() << "## **********************************************************";
	return;
}

void meta::print_info_rhmc(std::ostream* os, const Inputparameters& params)
{
	print_info_global(os, params);
	print_info_configs_io(os, params);
	print_info_prng_io(os, params);
	print_info_observables_rhmc_io(os, params);
	print_info_gauge(os, params);
	print_info_fermion(os, params);
	*os << "## **********************************************************" << endl;
	*os << "## RHMC parameters: "  << endl;
	*os << "##  " << endl;
	*os << "## Rational Approximations info:" << endl;
	*os << "##   - Generation of phi:" << endl;
	*os << "##       + x^(+" << params.get_num_tastes() << "/8)" << endl;
	*os << "##       + order = " << params.get_metro_approx_ord() << endl;
	*os << "##       + range = [" << params.get_approx_lower() << " , "
	    << params.get_approx_upper() << "]" << endl;
	*os << "##   - Molecular Dynamics:" << endl;
	*os << "##       + x^(-" << params.get_num_tastes() << "/4)" << endl;
	*os << "##       + order = " << params.get_md_approx_ord() << endl;
	*os << "##       + range = [" << params.get_approx_lower() << " , "
	    << params.get_approx_upper() << "]" << endl;
	*os << "##   - Evaluation of new action in Metropolis test:" << endl;
	*os << "##       + x^(-" << params.get_num_tastes() << "/4)" << endl;
	*os << "##       + order = " << params.get_metro_approx_ord() << endl;
	*os << "##       + range = [" << params.get_approx_lower() << " , "
	    << params.get_approx_upper() << "]" << endl;
	*os << "##  " << endl;
	*os << "## Strategy for finding max and min MdagM eigenvalue: " << (params.get_conservative() ? "conservative" : "NOT conservative") << endl;
	*os << "##  " << endl;
	*os << "## Simulation info:" << endl;
	*os << "##   - RHMC steps  = " << params.get_rhmcsteps() << endl;
	*os << "##   - precision used in  Mol. Dyn.-inversions = " << params.get_force_prec() << endl;
	*os << "##   - precision used in Metropolis-inversions = " << params.get_solver_prec() << endl;
	*os << "##  " << endl;
	*os << "## MD integration info:" << endl;
	*os << "##   - Time of integration  = " << params.get_tau() << endl;
	*os << "##   - Number of Timescales  = " << params.get_num_timescales() << endl;
	if(params.get_use_gauge_only() && params.get_num_timescales() == 1) {
		*os << "## use PureGaugeTheory in RHMC!" << endl;
	} else if (params.get_use_gauge_only() && params.get_num_timescales() > 1) {
		*os << "PureGaugeTheory can only be used with one timescale!\nPlease change the input-file!\nAborting..." << endl;
		throw out_of_range("Too many timescales");
	}
	//integrator infos
	*os << "##  " << endl;
	*os << "## Integrators info:" << endl;
	for(int i = 0; i < params.get_num_timescales(); i++) {
		print_info_integrator(os, i, params);
	}
	if(params.get_use_mp() == true) {
		//Not yet implemented!
		/*
		*os << "##  " << endl;
		*os <<  "## use mass preconditioning:" << endl;
		if(params.get_fermact_mp() == Inputparameters::wilson) {
		  *os <<  "## mp action: unimproved Wilson" << endl;
		  *os << "## kappa_mp  = " << params.get_kappa_mp() << endl;
		}
		if(params.get_fermact_mp() == Inputparameters::twistedmass) {
		  *os <<  "## mp action: twisted mass Wilson" << endl;
		  *os << "## kappa_mp  = " << params.get_kappa_mp() << endl;
		  *os << "## mu_mp     = " << params.get_mu_mp() << endl;
		}
		if(params.get_fermact_mp() == Inputparameters::clover) {
		  *os <<  "## mp action: clover Wilson" << endl;
		  *os << "## kappa_mp  = " << params.get_kappa_mp() << endl;
		  *os << "## csw_mp   = " << params.get_csw_mp() << endl;
		}
		*os << "##"  << endl;
		switch(params.get_solver_mp()) {
		  case Inputparameters::cg:
		    *os << "## Use CG-solver for mp inversions"  << endl;
		    break;
		  case Inputparameters::bicgstab:
		    *os << "## Use BiCGStab for mp inversions" << endl;
		    break;
		  case Inputparameters::bicgstab_save:
		    *os << "## Use BiCGStab-SAVE for mp inversions" << endl;
		    break;
		}
		*os << "## cgmax_mp  = " << params.get_cgmax_mp() << endl;
		*os << "## iter_refresh_mp  = " << params.get_iter_refresh_mp() << endl;
		*os << "##"  << endl;
		*/
	}
	*os << "## **********************************************************" << endl;
	return;
}

void meta::print_info_configs_io(const meta::Inputparameters& params)
{
	using namespace meta;

	logger.info() << "## **********************************************************";
	logger.info() << "## Configuration naming parameters:";
	logger.info() << "## digits in name:  " << params.get_config_number_digits();
	logger.info() << "## name prefix:   " << params.get_config_prefix();
	logger.info() << "## name postfix:   " << params.get_config_postfix();
	if(params.get_read_multiple_configs() == false) {
		logger.info() << "## read multiple configs: off";
	} else {
		logger.info() << "## read multiple configs: on";
		logger.info() << "## start nr:  " << params.get_config_read_start();
		logger.info() << "## end nr:   " << params.get_config_read_end();
		logger.info() << "## increment:    " << params.get_config_read_incr();
	}
}

void meta::print_info_configs_io(std::ostream * os, const meta::Inputparameters& params)
{
	using namespace meta;

	*os << "## **********************************************************" << endl;
	*os << "## Configuration naming parameters:" << endl;
	*os << "## digits in name:  " << params.get_config_number_digits() << endl;
	*os << "## name prefix:   " << params.get_config_prefix() << endl;
	*os << "## name postfix:   " << params.get_config_postfix() << endl;
	if(params.get_read_multiple_configs() == false) {
		*os << "## read multiple configs: off" << endl;
	} else {
		*os << "## read multiple configs: on" << endl;
		*os << "## start nr:  " << params.get_config_read_start() << endl;
		*os << "## end nr:   " << params.get_config_read_end() << endl;
		*os << "## increment:    " << params.get_config_read_incr() << endl;
	}
}

void meta::print_info_prng_io(const meta::Inputparameters& params)
{
	using namespace meta;

	logger.info() << "## **********************************************************";
	logger.info() << "## PRNG naming parameters:";
	logger.info() << "## digits in name:  " << params.get_config_number_digits();
	logger.info() << "## name prefix:   " << params.get_prng_prefix();
	logger.info() << "## name postfix:   " << params.get_prng_postfix();
}

void meta::print_info_prng_io(std::ostream * os, const meta::Inputparameters& params)
{
	using namespace meta;

	*os << "## **********************************************************" << endl;
	*os << "## Configuration naming parameters:" << endl;
	*os << "## digits in name:  " << params.get_config_number_digits() << endl;
	*os << "## name prefix:   " << params.get_prng_prefix() << endl;
	*os << "## name postfix:   " << params.get_prng_postfix() << endl;
}


void meta::print_info_observables_gauge_io(const meta::Inputparameters& params)
{
	using namespace meta;

	logger.info() << "## **********************************************************";
	logger.info() << "## gauge observables file naming parameters:";
	logger.info() << "## name prefix:   " << params.get_gauge_obs_prefix();
	logger.info() << "## name postfix:   " << params.get_gauge_obs_postfix();
	if(params.get_gauge_obs_to_single_file() == true) {
		logger.info() << "## write gauge observables to single file";
	} else {
		logger.info() << "## write gauge observables to multiple files";
	}
}

void meta::print_info_observables_gauge_io(std::ostream * os, const meta::Inputparameters& params)
{
	using namespace meta;

	*os << "## **********************************************************" << endl;
	*os << "## gauge observables file naming parameters:" << endl;
	*os << "## name prefix:   " << params.get_gauge_obs_prefix() << endl;
	*os << "## name postfix:   " << params.get_gauge_obs_postfix() << endl;
	if(params.get_gauge_obs_to_single_file() == true) {
		*os << "## write gauge observables to single file" << endl;
	} else {
		*os << "## write gauge observables to multiple files" << endl;
	}
}

static void print_info_observables_fermion_io(const meta::Inputparameters& params)
{
	using namespace meta;

	logger.info() << "## **********************************************************";
	logger.info() << "## fermionic observables file naming parameters:";
	if (params.get_measure_correlators() == true) {
		logger.info() << "## measure correlators in direction: " << params.get_corr_dir();
		logger.info() << "## correlators name prefix:   " << params.get_ferm_obs_corr_prefix();
		logger.info() << "## correlators name postfix:   " << params.get_ferm_obs_corr_postfix();
	}
	if (params.get_measure_pbp() == true) {
		logger.info() << "## chiral condensate name prefix:   " << params.get_ferm_obs_pbp_prefix();
		logger.info() << "## chiral condensate name postfix:   " << params.get_ferm_obs_pbp_postfix();
		if(params.get_pbp_version() == common::pbp_version::std )
			logger.info() << "## measure chiral condensate in standard version";
		if(params.get_pbp_version() == common::pbp_version::tm_one_end_trick ) {
			logger.info() << "## measure chiral condensate in twisted-mass one end trick version";
			if(params.get_fermact() != common::action::twistedmass)
				logger.fatal() << "## using the one end trick without twisted-mass action!";
		}
		if(params.get_sourcetype() == common::point)
			logger.warn() << "## calculating chiral condensate without stochastic estimators!";
	}
	if (params.get_measure_pbp() == false && (params.get_measure_correlators() == false )) {
		logger.info() << "## do not measure fermionic observables!";
	}
	if(params.get_ferm_obs_to_single_file() == true) {
		logger.info() << "## write fermion observables to single file";
	} else {
		logger.info() << "## write fermion observables to multiple files";
	}
}

static void print_info_observables_fermion_io(std::ostream * os, const meta::Inputparameters& params)
{
	using namespace meta;

	*os << "## **********************************************************" << endl;
	*os << "## fermionic observables file naming parameters:" << endl;
	if (params.get_measure_correlators() == true) {
		*os << "## measure correlators in direction: " << params.get_corr_dir() << endl;
		*os << "## correlators name prefix:   " << params.get_ferm_obs_corr_prefix() << endl;
		*os << "## correlators name postfix:   " << params.get_ferm_obs_corr_postfix() << endl;
	}
	if (params.get_measure_pbp() == true) {
		*os << "## chiral condensate name prefix:   " << params.get_ferm_obs_pbp_prefix() << endl;
		*os << "## chiral condensate name postfix:   " << params.get_ferm_obs_pbp_postfix() << endl;
		if(params.get_pbp_version() == common::pbp_version::std )
			*os << "## measure chiral condensate in standard version" << endl;
		if(params.get_pbp_version() == common::pbp_version::tm_one_end_trick ) {
			*os << "## measure chiral condensate in twisted-mass one end trick version" << endl;
			if(params.get_fermact() != common::action::twistedmass)
				*os << "## using the one end trick without twisted-mass action!" << endl;
		}
		if(params.get_sourcetype() == common::point)
			*os << "## calculating chiral condensate without stochastic estimators!" << endl;
	}
	if (params.get_measure_pbp() == false && (params.get_measure_correlators() == false )) {
		*os << "## do not measure fermionic observables!" << endl;
	}
	if(params.get_ferm_obs_to_single_file() == true) {
		*os << "## write fermionic observables to single file" << endl;
	} else {
		*os << "## write fermionic observables to multiple files" << endl;
	}
}

static void print_info_observables_hmc_io(const meta::Inputparameters& params)
{
	using namespace meta;

	logger.info() << "## **********************************************************";
	logger.info() << "## hmc observables file naming parameters:";
	logger.info() << "## name prefix:   " << params.get_hmc_obs_prefix();
	logger.info() << "## name postfix:   " << params.get_hmc_obs_postfix();
	if(params.get_hmc_obs_to_single_file() == true) {
		logger.info() << "## write hmc observables to single file";
	} else {
		logger.info() << "## write hmc observables to multiple files";
	}
}

static void print_info_observables_hmc_io(std::ostream * os, const meta::Inputparameters& params)
{
	using namespace meta;

	*os << "## **********************************************************" << endl;
	*os << "## hmc observables file naming parameters:" << endl;
	*os << "## name prefix:   " << params.get_hmc_obs_prefix() << endl;
	*os << "## name postfix:   " << params.get_hmc_obs_postfix() << endl;
	if(params.get_hmc_obs_to_single_file() == true) {
		*os << "## write hmc observables to single file" << endl;
	} else {
		*os << "## write hmc observables to multiple files" << endl;
	}
}

static void print_info_observables_rhmc_io(const meta::Inputparameters& params)
{
	using namespace meta;

	logger.info() << "## **********************************************************";
	logger.info() << "## RHMC observables file naming parameters:";
	logger.info() << "## name prefix:   " << params.get_rhmc_obs_prefix();
	logger.info() << "## name postfix:   " << params.get_rhmc_obs_postfix();
	if(params.get_rhmc_obs_to_single_file() == true) {
		logger.info() << "## write rhmc observables to single file";
	} else {
		logger.info() << "## write rhmc observables to multiple files";
	}
}

static void print_info_observables_rhmc_io(std::ostream * os, const meta::Inputparameters& params)
{
	using namespace meta;

	*os << "## **********************************************************" << endl;
	*os << "## RHMC observables file naming parameters:" << endl;
	*os << "## name prefix:   " << params.get_rhmc_obs_prefix() << endl;
	*os << "## name postfix:   " << params.get_rhmc_obs_postfix() << endl;
	if(params.get_rhmc_obs_to_single_file() == true) {
		*os << "## write rhmc observables to single file" << endl;
	} else {
		*os << "## write rhmc observables to multiple files" << endl;
	}
}

static void print_info_source(const meta::Inputparameters& params)
{
	logger.info() << "## **********************************************************";
	logger.info() << "## Source parameters:";
	logger.info() << "##";
	if(params.get_sourcetype() == common::sourcetypes::point) {
		logger.info() << "## Use pointsource for inversion" ;
		logger.info() << "## Position (x,y,z,t): " << params.get_source_x() << " " <<  params.get_source_y() << " " <<  params.get_source_z() << " " <<  params.get_source_t();
	} else if(params.get_sourcetype() == common::sourcetypes::volume) {
		logger.info() << "## Use volume sources for inversion" ;
		logger.info() << "## Number of sources: " << params.get_num_sources();
	} else if(params.get_sourcetype() == common::sourcetypes::timeslice) {
		logger.info() << "## Use timeslice sources for inversion" ;
		logger.info() << "## Use timeslice: " << params.get_source_t();
		logger.info() << "## Number of sources: " << params.get_num_sources();
	} else if(params.get_sourcetype() == common::sourcetypes::zslice) {
		logger.info() << "## Use zslice sources for inversion" ;
		logger.info() << "## Use zslice: " << params.get_source_z();
		logger.info() << "## Number of sources: " << params.get_num_sources();
	}
	if(params.get_sourcecontent() == common::sourcecontents::one) {
		logger.info() << "## fill sources with one";
	}  else if(params.get_sourcecontent() == common::sourcecontents::z4) {
		logger.info() << "## fill sources with z4 noise";
	} else if(params.get_sourcecontent() == common::sourcecontents::gaussian) {
		logger.info() << "## fill sources with gaussian noise";
	}
	if(params.get_measure_correlators() && params.get_fermact() == common::action::rooted_stagg){
	    if(params.get_sourcetype() == common::sourcetypes::point && params.get_num_sources() != 3) {
	        logger.fatal() << "## Number of sources: " << params.get_num_sources();
	        throw Print_Error_Message("Pointsource with number of sources different than \"3\" is chosen. This is will not give the full point-to-all propagator!");
	    }
	}
    if(params.get_measure_correlators() && params.get_fermact() != common::action::rooted_stagg){
        if(params.get_sourcetype() == common::sourcetypes::point && params.get_num_sources() != 12) {
            logger.fatal() << "## Number of sources: " << params.get_num_sources();
            throw Print_Error_Message("Pointsource with number of sources different than \"12\" is chosen. This is will not give the full point-to-all propagator!");
        }
    }
	if(params.get_sourcetype() == common::sourcetypes::point && params.get_sourcecontent() != common::sourcecontents::one) {
		logger.warn() << "## Pointsource with content different than \"one\" is chosen. This is not implemented yet and has no effect!";
	}
}

static void print_info_source(std::ostream * os, const meta::Inputparameters& params)
{
	*os << "## **********************************************************" << endl;
	*os << "## Source parameters:" << endl;
	*os << "##" << endl;
	if(params.get_sourcetype() == common::sourcetypes::point) {
		*os << "## Use pointsource for inversion" << endl;
		*os << "## Position (x,y,z,t): " << params.get_source_x() << " " <<  params.get_source_y() << " " <<  params.get_source_z() << " " <<  params.get_source_t() << endl;
	} else if(params.get_sourcetype() == common::sourcetypes::volume) {
		*os << "## Use volume sources for inversion" << endl;
		*os << "## Number of sources: " << params.get_num_sources() << endl;
	} else if(params.get_sourcetype() == common::sourcetypes::timeslice) {
		*os << "## Use timeslice sources for inversion" << endl;
		*os << "## Use timeslice: " << params.get_source_t() << endl;
		*os << "## Number of sources: " << params.get_num_sources() << endl;
	} else if(params.get_sourcetype() == common::sourcetypes::zslice) {
		*os << "## Use zslice sources for inversion" << endl;
		*os << "## Use zslice: " << params.get_source_z() << endl;
		*os << "## Number of sources: " << params.get_num_sources() << endl;
	}
	if(params.get_sourcecontent() == common::sourcecontents::one) {
		*os << "## fill sources with one" << endl;
	}  else if(params.get_sourcecontent() == common::sourcecontents::z4) {
		*os << "## fill sources with z4 noise" << endl;
	} else if(params.get_sourcecontent() == common::sourcecontents::gaussian) {
		*os << "## fill sources with gaussian noise" << endl;
	}
	if(params.get_sourcetype() == common::sourcetypes::point && params.get_num_sources() != 12) {
		*os << "## Pointsource with number of sources different than \"12\" is chosen. This is will not give the full point-to-all propagator!" << endl;
		*os << "## Number of sources: " << params.get_num_sources() << endl;
	}
	if(params.get_sourcetype() == common::sourcetypes::point && params.get_sourcecontent() != common::sourcecontents::one) {
		*os << "## Pointsource with content different than \"one\" is chosen. This is not implemented yet and has no effect!" << endl;
	}
}

void meta::print_info_flavour_doublet_correlators(const meta::Inputparameters& params)
{
	using namespace meta;

	logger.info() << "# flavour doublet correlators";
	if(params.get_corr_dir() == 3) {
	  logger.info() << "# format: J P (for each interpolator (if different) ) z real complex"  ;
		logger.info() << "# (J = Spin (0 or 1), P = Parity (0 positive, 1 negative), z spatial distance, value (aggregate x y z)" ;
	} else {
	  logger.info() << "# format: J P (for each interpolator (if different) ) t real complex"  ;
		logger.info() << "# (J = Spin (0 or 1), P = Parity (0 positive, 1 negative), t timelike distance, value (aggregate x y z)" ;
	}
}

void meta::print_info_flavour_doublet_correlators(std::ostream * os, const meta::Inputparameters& params)
{
	using namespace meta;

	*os << "# flavour doublet correlators" << std::endl;
	if(params.get_corr_dir() == 3) {
		*os << "# format: J P z real complex"  << std::endl;
		*os << "# (J = Spin (0 or 1), P = Parity (0 positive, 1 negative), z spatial distance, value (aggregate x y z)" << std::endl;
	} else {
		*os << "# format: J P t real complex"  << std::endl;
		*os << "# (J = Spin (0 or 1), P = Parity (0 positive, 1 negative), t timelike distance, value (aggregate x y z)" << std::endl;
	}
}

std::string meta::createLogfileName(const char* name)
{
	return std::string(name) + std::string(".log");
}

