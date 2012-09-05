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
	if (params.get_startcondition() == START_FROM_SOURCE) {
		std::string sf = params.get_sourcefile();
		logger.info() << "## sourcefile = " << sf;
	}
	if (params.get_startcondition() == COLD_START) {
		logger.info() << "## cold start";
	}
	if (params.get_startcondition() == HOT_START) {
		logger.info() << "## hot start";
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
		*os << "## USE ANISOTROPY: ON";
	} else {
		*os << "## USE ANISOTROPY: OFF";
	}
	*os  << "## Number of devices demanded for calculations: " << params.get_device_count()  << endl;
	*os  << "## **********************************************************" << endl;
	*os  << "## I/O parameters:" << endl;
	*os  << "## SvConf:  " << params.get_saveconfigs() << endl;
	*os  << "## WrFreq:  " << params.get_writefrequency() << endl;
	*os  << "## SvFreq:  " << params.get_savefrequency() << endl;
	if (params.get_startcondition() == START_FROM_SOURCE) {
		string sf = params.get_sourcefile();
		*os  << "## sourcefile = " << sf << endl;
	}
	if (params.get_startcondition() == COLD_START) {
		*os  << "## cold start" << endl;
	}
	if (params.get_startcondition() == HOT_START) {
		*os  << "## hot start" << endl;
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

void meta::print_info_heatbath(const char* progname, std::ostream* os, const Inputparameters& params)
{
	*os  << "## Starting heatbath program, executable name: " << progname << endl;
	print_info_global(os, params);
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
	if(params.get_use_pointsource() == true) {
		logger.info() << "## Use pointsource for inversion" ;
		logger.info() << "## Position (x,y,z,t): " << params.get_pointsource_x() << " " <<  params.get_pointsource_y() << " " <<  params.get_pointsource_z() << " " <<  params.get_pointsource_t();
	}
	if(params.get_use_pointsource() == false) {
		logger.info() << "## Use stochastic sources for inversion" ;
		logger.info() << "## Number of sources: " << params.get_num_sources();
	}
	if(params.get_use_cg() == true)
		logger.info() << "## Use CG-solver for inversions" ;
	if(params.get_use_cg() == false) {
		if(params.get_use_bicgstab_save() == false)
			logger.info() << "## Use BiCGStab for inversions";
		else
			logger.info() << "## Use BiCGStab-SAVE for inversions";
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
	if(params.get_use_pointsource() == true) {
		*os  << "## Use pointsource for inversion"  << endl;
		logger.info() << "## Position (x,y,z,t): " << params.get_pointsource_x() << " " <<  params.get_pointsource_y() << " " <<  params.get_pointsource_z() << " " <<  params.get_pointsource_t();
	}
	if(params.get_use_pointsource() == false) {
		*os  << "## Use stochastic sources for inversion"  << endl;
		*os  << "## Number of sources: " << params.get_num_sources()  << endl;
	}
	if(params.get_use_cg() == true)
		*os << "## Use CG-solver for inversions"  << endl;
	if(params.get_use_cg() == false) {
		if(params.get_use_bicgstab_save() == false)
			*os << "## Use BiCGStab for inversions" << endl;
		else
			*os << "## Use BiCGStab-SAVE for inversions" << endl;
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
	print_info_fermion(params);
	logger.info() << "## **********************************************************";
	return;
}

void meta::print_info_inverter(const char* progname, std::ostream* os, const Inputparameters& params)
{
	*os << "## Starting inverter program, executable name: " << progname << endl;
	print_info_global(os, params);
	print_info_fermion(os, params);
	*os << "## **********************************************************" << endl;
	return;
}

static void print_info_integrator(int number, const meta::Inputparameters& params)
{
	using namespace meta;

	string integrator_name;
	bool print_lambda = false;
	if(params.get_integrator(number) == LEAPFROG)
		integrator_name = "LEAPFROG";
	else if (params.get_integrator(number) == TWOMN) {
		integrator_name = "2MN";
		print_lambda = true;
	} else {
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
	if(params.get_integrator(number) == LEAPFROG)
		integrator_name = "LEAPFROG";
	else if (params.get_integrator(number) == TWOMN) {
		integrator_name = "2MN";
		print_lambda = true;
	} else {
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
		if(params.get_use_cg_mp() == true)
			logger.info() << "## Use CG-solver for mp inversions" ;
		if(params.get_use_cg_mp() == false) {
			if(params.get_use_bicgstab_save_mp() == false)
				logger.info() << "## Use BiCGStab for mp inversions";
			else
				logger.info() << "## Use BiCGStab-SAVE for mp inversions";
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
		if(params.get_use_cg_mp() == true)
			*os << "## Use CG-solver for mp inversions"  << '\n';
		if(params.get_use_cg_mp() == false) {
			if(params.get_use_bicgstab_save_mp() == false)
				*os << "## Use BiCGStab for mp inversions"  << '\n';
			else
				*os << "## Use BiCGStab-SAVE for mp inversions"  << '\n';
		}
		*os << "## cgmax_mp  = " << params.get_cgmax_mp()  << '\n';
		*os << "## iter_refresh_mp  = " << params.get_iter_refresh_mp()  << '\n';
		*os << "##"  << '\n';
	}
	*os << "## **********************************************************" << '\n';
	return;
}

