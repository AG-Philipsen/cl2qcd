/** @file
 * IO utility functions
 */

#include <stdexcept>
#include "gitcommitid.h"
#include "util.hpp"
#include "../logger.hpp"

using namespace std;

static void print_info_global(const meta::Inputparameters& params);
static void print_info_global(std::ostream* os, const meta::Inputparameters& params);
static void print_info_fermion(const meta::Inputparameters& params);
static void print_info_fermion(std::ostream * os, const meta::Inputparameters& params);
static void print_info_gauge(std::ostream* os, const meta::Inputparameters& params);
static void print_info_gauge(const meta::Inputparameters& params);
static void print_info_integrator(int number, const meta::Inputparameters& params);
static void print_info_integrator(std::ostream* os, int number, const meta::Inputparameters& params);
static void print_info_configs_io(const meta::Inputparameters& params);
static void print_info_configs_io(std::ostream * os, const meta::Inputparameters& params);
static void print_info_observables_gauge_io(const meta::Inputparameters& params);
static void print_info_observables_gauge_io(std::ostream * os, const meta::Inputparameters& params);
static void print_info_observables_fermion_io(const meta::Inputparameters& params);
static void print_info_observables_fermion_io(std::ostream * os, const meta::Inputparameters& params);
static void print_info_observables_hmc_io(const meta::Inputparameters& params);
static void print_info_observables_hmc_io(std::ostream * os, const meta::Inputparameters& params);
static void print_info_source(const meta::Inputparameters params);
static void print_info_source(std::ostream * os, const meta::Inputparameters params);

static void print_info_global(const meta::Inputparameters& params)
{
	using namespace meta;

	logger.info() << "## Build based on commit: " << GIT_COMMIT_ID;
	logger.info() << "## **********************************************************";
	logger.info() << "## Global parameters:";
	logger.info() << "## NSPACE:  " << params.get_nspace();
	logger.info() << "## NTIME:   " << params.get_ntime();
	logger.info() << "## NDIM:    " << NDIM;
	logger.info() << "## NCOLOR:  " << NC;
	logger.info() << "## NSPIN:   " << NSPIN;
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
	logger.info() << "## **********************************************************";
	logger.info() << "## I/O parameters:";
	logger.info() << "## SvConf:  " << params.get_saveconfigs();
	logger.info() << "## WrFreq:  " << params.get_writefrequency();
	logger.info() << "## SvFreq:  " << params.get_savefrequency();
	switch(params.get_startcondition()) {
		case Inputparameters::start_from_source: {
			std::string sf = params.get_sourcefile();
			logger.info() << "## sourcefile = " << sf;
		}
		break;
		case Inputparameters::cold_start:
			logger.info() << "## cold start";
			break;
		case Inputparameters::hot_start:
			logger.info() << "## hot start";
			break;
	}
	if(params.get_use_smearing() == 1) {
		logger.info() << "## **********************************************************";
		logger.info() << "## Apply Smearing with:";
		logger.info() << "## rho:      " << params.get_rho();
		logger.info() << "## rho_iter: " << params.get_rho_iter();
	}
}

static void print_info_global(std::ostream* os, const meta::Inputparameters& params)
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
	*os  << "## **********************************************************" << endl;
	*os  << "## I/O parameters:" << endl;
	*os  << "## SvConf:  " << params.get_saveconfigs() << endl;
	*os  << "## WrFreq:  " << params.get_writefrequency() << endl;
	*os  << "## SvFreq:  " << params.get_savefrequency() << endl;
	switch(params.get_startcondition()) {
		case Inputparameters::start_from_source: {
			std::string sf = params.get_sourcefile();
			*os << "## sourcefile = " << sf << endl;;
		}
		break;
		case Inputparameters::cold_start:
		  *os << "## cold start" << endl;;
			break;
		case Inputparameters::hot_start:
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


void meta::print_info_heatbath(const char* progname, const Inputparameters& params)
{
	logger.info() << "## Starting heatbath program, executable name: " << progname;
	print_info_global(params);
	print_info_configs_io(params);
	print_info_observables_gauge_io(params);
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

void meta::print_info_gaugeobservables(const char* progname, std::ostream* os, const Inputparameters& params)
{
	*os  << "## Starting gaugeobservables program, executable name: " << progname << endl;
	print_info_global(os, params);
	print_info_configs_io(os, params);
	return;
}

void meta::print_info_gaugeobservables(const char* progname, const Inputparameters& params)
{
	logger.info() << "## Starting gaugeobservables program, executable name: " << progname;
	print_info_global(params);
	print_info_configs_io(params);
	return;
}

void meta::print_info_heatbath(const char* progname, std::ostream* os, const Inputparameters& params)
{
	*os  << "## Starting heatbath program, executable name: " << progname << endl;
	print_info_global(os, params);
	print_info_configs_io(os, params);
	print_info_observables_gauge_io(os, params);
	*os  << "## **********************************************************" << endl;
	*os  << "## Simulation parameters:" << endl;
	*os  << "## beta           = " << params.get_beta() << endl;
	*os  << "## xi             = " << params.get_xi();
	*os  << "## thermsteps     = " << params.get_thermalizationsteps() << endl;
	*os  << "## heatbathsteps  = " << params.get_heatbathsteps() << endl;
	*os  << "## overrelaxsteps = " << params.get_overrelaxsteps() << endl;
	*os  << "## **********************************************************" << endl;
	return;
}

void meta::print_info_tkkappa(const char* progname, std::ostream* os, const Inputparameters& params)
{
	*os << "## Starting tk_kappa program, " << progname << endl;
	print_info_global(os, params);
	print_info_configs_io(os, params);
	print_info_observables_gauge_io(os, params);
	*os  << "## **********************************************************" << endl;
	*os  << "## Simulation parameters:" << endl;
	*os  << "## beta           = " << params.get_beta() << endl;
	*os  << "## xi             = " << params.get_xi();
	*os  << "## thermsteps     = " << params.get_thermalizationsteps() << endl;
	*os  << "## heatbathsteps  = " << params.get_heatbathsteps() << endl;
	*os  << "## overrelaxsteps = " << params.get_overrelaxsteps() << endl;
	*os  << "## TODO: INSERT SPECIFIC PARAMETERS!!!!!" << endl;
	*os  << "## **********************************************************" << endl;
	return;
}

void meta::print_info_tkkappa(const char* progname, const Inputparameters& params)
{
	logger.info() << "## Starting tk_kappa program, " << progname ;
	print_info_global(params);
	print_info_configs_io(params);
	print_info_observables_gauge_io(params);
	logger.info() << "## **********************************************************";
	logger.info() << "## Simulation parameters:";
	logger.info() << "## beta           = " << params.get_beta();
	logger.info() << "## xi             = " << params.get_xi();
	logger.info() << "## thermsteps     = " << params.get_thermalizationsteps() ;
	logger.info() << "## heatbathsteps  = " << params.get_heatbathsteps();
	logger.info() << "## overrelaxsteps = " << params.get_overrelaxsteps();
	logger.info() << "## TODO: INSERT SPECIFIC PARAMETERS!!!!!";
	logger.info() << "## **********************************************************";
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
	if(params.get_fermact() == Inputparameters::wilson) {
		logger.info() <<  "## fermion action: unimproved Wilson";
		logger.info() << "## kappa  = " << params.get_kappa();
	}
	if(params.get_fermact() == Inputparameters::twistedmass) {
		logger.info() <<  "## fermion action: twisted mass Wilson";
		logger.info() << "## kappa  = " << params.get_kappa();
		logger.info() << "## mu     = " << params.get_mu();
	}
	if(params.get_fermact() == Inputparameters::clover) {
		logger.info() <<  "## fermion action: clover Wilson";
		logger.info() << "## kappa  = " << params.get_kappa();
		logger.info() << "## csw    = " << params.get_csw();
	}
	logger.info() << "##" ;
	logger.info() << "## Inverter parameters:";
	logger.info() << "## precision for inversions = " << params.get_solver_prec();
	if(params.get_use_eo() == true)
		logger.info() << "## Use even-odd preconditioning" ;
	if(params.get_use_eo() == false)
		logger.info() << "## Do NOT use even-odd preconditioning";
	switch(params.get_solver()) {
		case Inputparameters::cg:
			logger.info() << "## Use CG-solver for inversions" ;
			break;
		case Inputparameters::bicgstab:
			logger.info() << "## Use BiCGStab for inversions";
			break;
		case Inputparameters::bicgstab_save:
			logger.info() << "## Use BiCGStab-SAVE for inversions";
			break;
	}
	logger.info() << "## cgmax  = " << params.get_cgmax();
	logger.info() << "## iter_refresh  = " << params.get_iter_refresh();

	if(params.get_profile_solver() == true)
		logger.warn() << "## Profiling of solver activated. This may influence the overall performance time!";

	if(params.get_use_merge_kernels_fermion() == true)
		logger.info() << "## Use merged fermionmatrix kernels where implemented!!";

	if(params.get_use_merge_kernels_spinor() == true)
		logger.info() << "## Use merged spinor kernels where implemented!!";

	//print extra warning if BC are set to default since this is a serious source of errors...
	if ( params.get_theta_fermion_spatial() == 0. && params.get_theta_fermion_temporal() == 0.) {
		logger.warn() << "\nNOTE: BCs have been set to periodic values by default!!\nTo change this use e.g. ThetaT/ThetaS in the input-file.\n";
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
	if(params.get_fermact() == Inputparameters::wilson) {
		*os <<  "## fermion action: unimproved Wilson" << endl;
		*os  << "## kappa  = " << params.get_kappa() << endl;
	}
	if(params.get_fermact() == Inputparameters::twistedmass) {
		*os <<  "## fermion action: twisted mass Wilson" << endl;
		*os  << "## kappa  = " << params.get_kappa() << endl;
		*os  << "## mu     = " << params.get_mu() << endl;
	}
	if(params.get_fermact() == Inputparameters::clover) {
		*os <<  "## fermion action: clover Wilson" << endl;
		*os  << "## kappa  = " << params.get_kappa() << endl;
		*os  << "## csw    = " << params.get_csw() << endl;
	}
	*os  << "##" << endl;
	*os  << "## Inverter parameters:" << endl;
	*os << "## precision for inversions = " << params.get_solver_prec() << endl;
	if(params.get_use_eo() == true)
		*os  << "## Use even-odd preconditioning" << endl;
	if(params.get_use_eo() == false)
		*os  << "## Do NOT use even-odd preconditioning" << endl;
	switch(params.get_solver()) {
		case Inputparameters::cg:
			*os << "## Use CG-solver for inversions" ;
			break;
		case Inputparameters::bicgstab:
			*os << "## Use BiCGStab for inversions";
			break;
		case Inputparameters::bicgstab_save:
			*os << "## Use BiCGStab-SAVE for inversions";
			break;
	}
	*os << "## cgmax  = " << params.get_cgmax() << endl;
	*os << "## iter_refresh  = " << params.get_iter_refresh() << endl;

	if(params.get_profile_solver() == true)
		*os << "## Profiling of solver activated. This may influence the overall performance time!" << endl;

	if(params.get_use_merge_kernels_fermion() == true)
		*os << "## Use merged fermionmatrix kernels where implemented!!" << endl;

	if(params.get_use_merge_kernels_spinor() == true)
		*os << "## Use merged spinor kernels where implemented!!" << endl;


	//print extra warning if BC are set to default since this is a serious source of errors...
	if ( params.get_theta_fermion_spatial() == 0. && params.get_theta_fermion_temporal() == 0.) {
		*os << "\nNOTE: BCs have been set to periodic values by default!!\nTo change this use e.g. ThetaT/ThetaS in the input-file.\n";
	}
}

static void print_info_gauge(std::ostream* os, const meta::Inputparameters& params)
{
	using namespace meta;

	*os << "## **********************************************************" << endl;
	*os << "## Gauge parameters:" << endl;
	*os << "##" << endl;
	if(params.get_gaugeact() == Inputparameters::wilson) {
		*os <<  "## gauge action: unimproved Wilson" << endl;
	}
	if(params.get_gaugeact() == Inputparameters::tlsym) {
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
	if(params.get_gaugeact() == Inputparameters::wilson) {
		logger.info() <<  "## gauge action: unimproved Wilson";
	}
	if(params.get_gaugeact() == Inputparameters::tlsym) {
		logger.info() <<  "## gauge action: tree level Symanzik";
		logger.info() << "## c0  = " << get_c0(params);
		logger.info() << "## c1  = " << get_c1(params);
	}
}

void meta::print_info_inverter(const char* progname, const Inputparameters& params)
{
	logger.info() << "## Starting inverter program, executable name: " << progname;
	print_info_global(params);
	print_info_configs_io(params);
	print_info_observables_fermion_io(params);
	print_info_fermion(params);
	print_info_source(params);
	logger.info() << "## **********************************************************";
	return;
}

void meta::print_info_inverter(const char* progname, std::ostream* os, const Inputparameters& params)
{
	*os << "## Starting inverter program, executable name: " << progname << endl;
	print_info_global(os, params);
	print_info_configs_io(os, params);
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
		case Inputparameters::leapfrog:
			integrator_name = "LEAPFROG";
			break;
		case Inputparameters::twomn:
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
		case Inputparameters::leapfrog:
			integrator_name = "LEAPFROG";
			break;
		case Inputparameters::twomn:
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

void meta::print_info_hmc(const char* progname, const Inputparameters& params)
{

	logger.info() << "## Starting hmc program, executable name: " << progname ;
	print_info_global(params);
	print_info_configs_io(params);
	print_info_observables_hmc_io(params);
	print_info_fermion(params);
	print_info_gauge(params);
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
	}
	logger.info() << "## **********************************************************";
	return;
}

void meta::print_info_hmc(const char* progname, std::ostream* os, const Inputparameters& params)
{
	*os << "## Starting hmc program, executable name: " << progname << endl;
	print_info_global(os, params);
	print_info_configs_io(os, params);
	print_info_observables_hmc_io(os, params);
	print_info_fermion(os, params);
	print_info_gauge(os, params);
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
		if(params.get_fermact_mp() == Inputparameters::wilson) {
			*os <<  "## mp action: unimproved Wilson"  << '\n';
			*os << "## kappa_mp  = " << params.get_kappa_mp()  << '\n';
		}
		if(params.get_fermact_mp() == Inputparameters::twistedmass) {
			*os <<  "## mp action: twisted mass Wilson"  << '\n';
			*os << "## kappa_mp  = " << params.get_kappa_mp()  << '\n';
			*os << "## mu_mp     = " << params.get_mu_mp()  << '\n';
		}
		if(params.get_fermact_mp() == Inputparameters::clover) {
			*os <<  "## mp action: clover Wilson";
			*os << "## kappa_mp  = " << params.get_kappa_mp()  << '\n';
			*os << "## csw_mp   = " << params.get_csw_mp()  << '\n';
		}
		*os << "##"  << endl;
		switch(params.get_solver_mp()) {
			case Inputparameters::cg:
				*os << "## Use CG-solver for mp inversions" ;
				break;
			case Inputparameters::bicgstab:
				*os << "## Use BiCGStab for mp inversions";
				break;
			case Inputparameters::bicgstab_save:
				*os << "## Use BiCGStab-SAVE for mp inversions";
				break;
		}
		*os << "## cgmax_mp  = " << params.get_cgmax_mp()  << '\n';
		*os << "## iter_refresh_mp  = " << params.get_iter_refresh_mp()  << '\n';
		*os << "##"  << '\n';
	}
	*os << "## **********************************************************" << '\n';
	return;
}

static void print_info_configs_io(const meta::Inputparameters& params)
{
	using namespace meta;

	logger.info() << "## **********************************************************";
	logger.info() << "## CONFIGURATION NAMING PARAMETERS:";
	logger.info() << "## DIGITS IN NAME:  " << params.get_config_number_digits();
	logger.info() << "## NAME PREFIX:   " << params.get_config_prefix();
	logger.info() << "## NAME POSTFIX:   " << params.get_config_postfix();
	if(params.get_read_multiple_configs() == false) {
		logger.info() << "## READ MULTIPLE CONFIGS: OFF";
	} else {
		logger.info() << "## READ MULTIPLE CONFIGS: ON";
		logger.info() << "## START NR:  " << params.get_config_read_start();
		logger.info() << "## END NR:   " << params.get_config_read_end();
		logger.info() << "## INCREMENT:    " << params.get_config_read_incr();
	}
}

static void print_info_configs_io(std::ostream * os, const meta::Inputparameters& params)
{
	using namespace meta;

	*os << "## **********************************************************" << endl;
	*os << "## CONFIGURATION NAMING PARAMETERS:" << endl;
	*os << "## DIGITS IN NAME:  " << params.get_config_number_digits() << endl;
	*os << "## NAME PREFIX:   " << params.get_config_prefix() << endl;
	*os << "## NAME POSTFIX:   " << params.get_config_postfix() << endl;
	if(params.get_read_multiple_configs() == false) {
		*os << "## READ MULTIPLE CONFIGS: OFF" << endl;
	} else {
		*os << "## READ MULTIPLE CONFIGS: ON" << endl;
		*os << "## START NR:  " << params.get_config_read_start() << endl;
		*os << "## END NR:   " << params.get_config_read_end() << endl;
		*os << "## INCREMENT:    " << params.get_config_read_incr() << endl;
	}
}

static void print_info_observables_gauge_io(const meta::Inputparameters& params)
{
	using namespace meta;

	logger.info() << "## **********************************************************";
	logger.info() << "## GAUGE OBSERVABLES FILE NAMING PARAMETERS:";
	logger.info() << "## NAME PREFIX:   " << params.get_gauge_obs_prefix();
	logger.info() << "## NAME POSTFIX:   " << params.get_gauge_obs_postfix();
	if(params.get_gauge_obs_to_single_file() == true) {
		logger.info() << "## WRITE GAUGE OBSERVABLES TO SINGLE FILE";
	} else {
		logger.info() << "## WRITE GAUGE OBSERVABLES TO MULTIPLE FILES";
	}
}

static void print_info_observables_gauge_io(std::ostream * os, const meta::Inputparameters& params)
{
	using namespace meta;

	*os<< "## **********************************************************"<< endl;
	*os<< "## GAUGE OBSERVABLES FILE NAMING PARAMETERS:"<< endl;
	*os<< "## NAME PREFIX:   " << params.get_gauge_obs_prefix()<< endl;
	*os<< "## NAME POSTFIX:   " << params.get_gauge_obs_postfix()<< endl;
	if(params.get_gauge_obs_to_single_file() == true) {
		*os<< "## WRITE GAUGE OBSERVABLES TO SINGLE FILE"<< endl;
	} else {
		*os<< "## WRITE GAUGE OBSERVABLES TO MULTIPLE FILES"<< endl;
	}
}

static void print_info_observables_fermion_io(const meta::Inputparameters& params)
{
	using namespace meta;

	logger.info() << "## **********************************************************";
	logger.info() << "## FERMIONIC OBSERVABLES FILE NAMING PARAMETERS:";
	logger.info() << "## NAME PREFIX:   " << params.get_ferm_obs_prefix();
	logger.info() << "## NAME POSTFIX:   " << params.get_ferm_obs_postfix();
	if(params.get_ferm_obs_to_single_file() == true) {
		logger.info() << "## WRITE FERMION OBSERVABLES TO SINGLE FILE";
	} else {
		logger.info() << "## WRITE FERMION OBSERVABLES TO MULTIPLE FILES";
	}
}

static void print_info_observables_fermion_io(std::ostream * os, const meta::Inputparameters& params)
{
	using namespace meta;

	*os<< "## **********************************************************"<< endl;
	*os<< "## FERMIONIC OBSERVABLES FILE NAMING PARAMETERS:"<< endl;
	*os<< "## NAME PREFIX:   " << params.get_ferm_obs_prefix()<< endl;
	*os<< "## NAME POSTFIX:   " << params.get_ferm_obs_postfix()<< endl;
	if(params.get_ferm_obs_to_single_file() == true) {
		*os<< "## WRITE FERMIONIC OBSERVABLES TO SINGLE FILE"<< endl;
	} else {
		*os<< "## WRITE FERMIONIC OBSERVABLES TO MULTIPLE FILES"<< endl;
	}
}

static void print_info_observables_hmc_io(const meta::Inputparameters& params)
{
	using namespace meta;

	logger.info() << "## **********************************************************";
	logger.info() << "## HMC OBSERVABLES FILE NAMING PARAMETERS:";
	logger.info() << "## NAME PREFIX:   " << params.get_hmc_obs_prefix();
	logger.info() << "## NAME POSTFIX:   " << params.get_hmc_obs_postfix();
	if(params.get_hmc_obs_to_single_file() == true) {
		logger.info() << "## WRITE HMC OBSERVABLES TO SINGLE FILE";
	} else {
		logger.info() << "## WRITE HMC OBSERVABLES TO MULTIPLE FILES";
	}
}

static void print_info_observables_hmc_io(std::ostream * os, const meta::Inputparameters& params)
{
	using namespace meta;

	*os<< "## **********************************************************"<< endl;
	*os<< "## HMC OBSERVABLES FILE NAMING PARAMETERS:"<< endl;
	*os<< "## NAME PREFIX:   " << params.get_hmc_obs_prefix()<< endl;
	*os<< "## NAME POSTFIX:   " << params.get_hmc_obs_postfix()<< endl;
	if(params.get_hmc_obs_to_single_file() == true) {
		*os<< "## WRITE HMC OBSERVABLES TO SINGLE FILE"<< endl;
	} else {
		*os<< "## WRITE HMC OBSERVABLES TO MULTIPLE FILES"<< endl;
	}
}

static void print_info_source(const meta::Inputparameters params) 
{
    logger.info() << "## **********************************************************";
    logger.info() << "## Source parameters:";
    logger.info() << "##";
    if(params.get_sourcetype() == meta::Inputparameters::sourcetypes::point) {
	logger.info() << "## Use pointsource for inversion" ;
	logger.info() << "## Position (x,y,z,t): " << params.get_source_x() << " " <<  params.get_source_y() << " " <<  params.get_source_z() << " " <<  params.get_source_t();
    }
    else if(params.get_sourcetype() == meta::Inputparameters::sourcetypes::volume) {
      logger.info() << "## Use volume sources for inversion" ;
	logger.info() << "## Number of sources: " << params.get_num_sources();
	}
    else if(params.get_sourcetype() == meta::Inputparameters::sourcetypes::timeslice) {
	logger.info() << "## Use timeslice sources for inversion" ;
	logger.info() << "## Use timeslice: " << params.get_source_t();
	logger.info() << "## Number of sources: " << params.get_num_sources();
    }
    if(params.get_sourcecontent() == meta::Inputparameters::sourcecontents::one){
      logger.info() << "## fill sources with one";
    }  else if(params.get_sourcecontent() == meta::Inputparameters::sourcecontents::z4){
      logger.info() << "## fill sources with z4 noise";
    }
    else if(params.get_sourcecontent() == meta::Inputparameters::sourcecontents::gaussian){
      logger.info() << "## fill sources with gaussian noise";
    }
}

static void print_info_source(std::ostream * os, const meta::Inputparameters params) 
{
    *os<< "## **********************************************************"<< endl;
    *os<< "## Source parameters:"<< endl;
    *os<< "##" << endl;
    if(params.get_sourcetype() == meta::Inputparameters::sourcetypes::point) {
	*os<< "## Use pointsource for inversion" << endl;
	*os<< "## Position (x,y,z,t): " << params.get_source_x() << " " <<  params.get_source_y() << " " <<  params.get_source_z() << " " <<  params.get_source_t()<< endl;
    }
    else if(params.get_sourcetype() == meta::Inputparameters::sourcetypes::volume) {
      *os<< "## Use volume sources for inversion" << endl;
	*os<< "## Number of sources: " << params.get_num_sources()<< endl;
	}
    else if(params.get_sourcetype() == meta::Inputparameters::sourcetypes::timeslice) {
	*os<< "## Use timeslice sources for inversion" << endl;
	*os<< "## Number of sources: " << params.get_num_sources()<< endl;
	*os << "## Number of sources: " << params.get_num_sources() << endl;
    }
    if(params.get_sourcecontent() == meta::Inputparameters::sourcecontents::one){
      *os<< "## fill sources with one"<< endl;
    }  else if(params.get_sourcecontent() == meta::Inputparameters::sourcecontents::z4){
      *os<< "## fill sources with z4 noise"<< endl;
    }
    else if(params.get_sourcecontent() == meta::Inputparameters::sourcecontents::gaussian){
      *os<< "## fill sources with gaussian noise"<< endl;
    }
}

