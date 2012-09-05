/** @file
 * Testcases for the Inputparameters class
 *
 * (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 */

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE meta::Inputparameters
#include <boost/test/unit_test.hpp>

#include "inputparameters.hpp"

#include "../types.h"

using namespace meta;

void checkDefaults(const Inputparameters params)
{
	BOOST_REQUIRE_EQUAL(params.get_precision(), sizeof(double) * 8);

	BOOST_REQUIRE_EQUAL(params.get_selected_devices().size(), 0);
	BOOST_REQUIRE_EQUAL(params.get_device_count(), 0);

	BOOST_REQUIRE_EQUAL(params.get_use_aniso(), false);
	BOOST_REQUIRE_EQUAL(params.get_use_chem_pot_re(), false);
	BOOST_REQUIRE_EQUAL(params.get_use_chem_pot_im(), false);
	BOOST_REQUIRE_EQUAL(params.get_use_smearing(), false);
	BOOST_REQUIRE_EQUAL(params.get_use_mp(), false);
	BOOST_REQUIRE_EQUAL(params.get_nspace(), 4);
	BOOST_REQUIRE_EQUAL(params.get_ntime(), 8);

	BOOST_REQUIRE_EQUAL(params.get_startcondition(), Inputparameters::start_from_source);
	BOOST_REQUIRE_EQUAL(params.get_saveconfigs(), false);
	BOOST_REQUIRE_EQUAL(params.get_writefrequency(), 1);
	BOOST_REQUIRE_EQUAL(params.get_savefrequency(), 100);
	BOOST_REQUIRE_EQUAL(params.get_sourcefile(), "conf.00000");
	BOOST_REQUIRE_EQUAL(params.get_print_to_screen(), false);
	//This is obvious!!!
	BOOST_REQUIRE_EQUAL(params.get_host_seed(), 4815);

	//gaugefield parameters
	BOOST_REQUIRE_EQUAL(params.get_beta(), 4.0);
	BOOST_REQUIRE_EQUAL(params.get_rho(), 0.);
	BOOST_REQUIRE_EQUAL(params.get_rho_iter(), 0);
	BOOST_REQUIRE_EQUAL(params.get_gaugeact(), Inputparameters::wilson);

	//heatbath parameters
	BOOST_REQUIRE_EQUAL(params.get_thermalizationsteps(), 0);
	BOOST_REQUIRE_EQUAL(params.get_heatbathsteps(), 1000);
	BOOST_REQUIRE_EQUAL(params.get_overrelaxsteps(), 1);
	BOOST_REQUIRE_EQUAL(params.get_xi(), 1);

	//fermionic parameters
	BOOST_REQUIRE_EQUAL(params.get_fermact(), Inputparameters::wilson);
	BOOST_REQUIRE_EQUAL(params.get_fermact_mp(), Inputparameters::wilson);
	BOOST_REQUIRE_EQUAL(params.get_kappa(), 0.125);
	BOOST_REQUIRE_EQUAL(params.get_mu(), 0.006);
	BOOST_REQUIRE_EQUAL(params.get_csw(), 0.);
	BOOST_REQUIRE_EQUAL(params.get_iter0(), 0);
	BOOST_REQUIRE_EQUAL(params.get_iter1(), 0);
	BOOST_REQUIRE_EQUAL(params.get_kappa_mp(), 0.125);
	BOOST_REQUIRE_EQUAL(params.get_mu_mp(), 0.006);
	BOOST_REQUIRE_EQUAL(params.get_csw_mp(), 0.);
	BOOST_REQUIRE_EQUAL(params.get_iter0_mp(), 0);
	BOOST_REQUIRE_EQUAL(params.get_iter1_mp(), 0);
	BOOST_REQUIRE_EQUAL(params.get_cgmax(), 1000);
	BOOST_REQUIRE_EQUAL(params.get_cgmax_mp(), 1000);
	BOOST_REQUIRE_EQUAL(params.get_theta_fermion_spatial(), 0.);
	BOOST_REQUIRE_EQUAL(params.get_theta_fermion_temporal(), 0.);
	BOOST_REQUIRE_EQUAL(params.get_chem_pot_re(), 0.);
	BOOST_REQUIRE_EQUAL(params.get_chem_pot_im(), 0.);
	BOOST_REQUIRE_EQUAL(params.get_use_eo(), true);
	//at the moment, only 2 solvers are implemented..
	BOOST_REQUIRE_EQUAL(params.get_use_cg(), false);
	BOOST_REQUIRE_EQUAL(params.get_use_cg_mp(), false);
	BOOST_REQUIRE_EQUAL(params.get_use_bicgstab_save(), false);
	BOOST_REQUIRE_EQUAL(params.get_use_bicgstab_save_mp(), false);
	BOOST_REQUIRE_EQUAL(params.get_use_pointsource(), true);
	BOOST_REQUIRE_EQUAL(params.get_use_gauge_only(), false);
	BOOST_REQUIRE_EQUAL(params.get_num_sources(), 12);
	BOOST_REQUIRE_EQUAL(params.get_pointsource_x(), 0);
	BOOST_REQUIRE_EQUAL(params.get_pointsource_y(), 0);
	BOOST_REQUIRE_EQUAL(params.get_pointsource_z(), 0);
	BOOST_REQUIRE_EQUAL(params.get_pointsource_t(), 0);
#ifdef _USEDOUBLEPREC_
	BOOST_REQUIRE_EQUAL(params.get_solver_prec(), 1e-23);
	BOOST_REQUIRE_EQUAL(params.get_force_prec(), 1e-12);
#else
	BOOST_REQUIRE_EQUAL(params.get_solver_prec(), 1e-16);
	BOOST_REQUIRE_EQUAL(params.get_force_prec(), 1e-8);
#endif
	BOOST_REQUIRE_EQUAL(params.get_iter_refresh(), 100);
	BOOST_REQUIRE_EQUAL(params.get_iter_refresh_mp(), 100);

	//HMC specific parameters
	BOOST_REQUIRE_EQUAL(params.get_tau(), 0.5);
	BOOST_REQUIRE_EQUAL(params.get_reversibility_check(), false);
	BOOST_REQUIRE_EQUAL(params.get_integrationsteps(0), 10);
	BOOST_REQUIRE_EQUAL(params.get_integrationsteps(1), 10);
	BOOST_REQUIRE_EQUAL(params.get_integrationsteps(2), 10);
	BOOST_REQUIRE_EQUAL(params.get_hmcsteps(), 10);
	BOOST_REQUIRE_EQUAL(params.get_num_timescales(), 1);
	BOOST_REQUIRE_EQUAL(params.get_integrator(0), Inputparameters::leapfrog);
	BOOST_REQUIRE_EQUAL(params.get_integrator(1), Inputparameters::leapfrog);
	BOOST_REQUIRE_EQUAL(params.get_integrator(2), Inputparameters::leapfrog);
	//this is the optimal value...
	BOOST_REQUIRE_EQUAL(params.get_lambda(0), 0.1931833275037836);
	BOOST_REQUIRE_EQUAL(params.get_lambda(1), 0.1931833275037836);
	BOOST_REQUIRE_EQUAL(params.get_lambda(2), 0.1931833275037836);

	//direction for the correlator
	BOOST_REQUIRE_EQUAL(params.get_corr_dir(), 3);

	BOOST_REQUIRE_EQUAL(params.get_use_same_rnd_numbers(), false);
	BOOST_REQUIRE_EQUAL(params.get_profile_solver(), false);

	BOOST_REQUIRE_EQUAL(params.is_ocl_compiler_opt_disabled(), false);

	BOOST_REQUIRE_EQUAL(params.get_use_merge_kernels_spinor(), false);
	BOOST_REQUIRE_EQUAL(params.get_use_merge_kernels_fermion(), false);
}

BOOST_AUTO_TEST_CASE(defaults)
{
	const char* _params[] = {"foo"};
	Inputparameters params(1, _params);
	checkDefaults(params);
}

BOOST_AUTO_TEST_CASE(input_file1)
{
	// empty input file
	const char* _params[] = {"foo", "test_input_1"};
	Inputparameters params(2, _params);
	checkDefaults(params);
}

BOOST_AUTO_TEST_CASE(input_file2)
{
	// filled input file
	const char* _params2[] = {"foo", "test_input_2"};
	Inputparameters params(2, _params2);
	BOOST_REQUIRE_EQUAL(params.get_precision(), 32);
	BOOST_REQUIRE_EQUAL(params.get_selected_devices().size(), 0);
	BOOST_REQUIRE_EQUAL(params.get_device_count(), 1);

	BOOST_REQUIRE_EQUAL(params.get_use_aniso(), true);
	BOOST_REQUIRE_EQUAL(params.get_use_chem_pot_re(), true);
	BOOST_REQUIRE_EQUAL(params.get_use_chem_pot_im(), true);
	BOOST_REQUIRE_EQUAL(params.get_use_smearing(), true);
	BOOST_REQUIRE_EQUAL(params.get_use_mp(), true);
	BOOST_REQUIRE_EQUAL(params.get_nspace(), 32);
	BOOST_REQUIRE_EQUAL(params.get_ntime(), 12);

	BOOST_REQUIRE_EQUAL(params.get_startcondition(), Inputparameters::start_from_source);
	BOOST_REQUIRE_EQUAL(params.get_saveconfigs(), true);
	BOOST_REQUIRE_EQUAL(params.get_writefrequency(), 10);
	BOOST_REQUIRE_EQUAL(params.get_savefrequency(), 10);
	BOOST_REQUIRE_EQUAL(params.get_sourcefile(), "conf.11111");
	BOOST_REQUIRE_EQUAL(params.get_print_to_screen(), true);
	//This is obvious!!!
	BOOST_REQUIRE_EQUAL(params.get_host_seed(), 42);

	//gaugefield parameters
	BOOST_REQUIRE_EQUAL(params.get_beta(), 4.1);
	BOOST_REQUIRE_EQUAL(params.get_rho(), 0.1);
	BOOST_REQUIRE_EQUAL(params.get_rho_iter(), 4);
	BOOST_REQUIRE_EQUAL(params.get_gaugeact(), Inputparameters::twistedmass);

	//heatbath parameters
	BOOST_REQUIRE_EQUAL(params.get_thermalizationsteps(), 10);
	BOOST_REQUIRE_EQUAL(params.get_heatbathsteps(), 100);
	BOOST_REQUIRE_EQUAL(params.get_overrelaxsteps(), 10);
	BOOST_REQUIRE_EQUAL(params.get_xi(), 2);

	//fermionic parameters
	BOOST_REQUIRE_EQUAL(params.get_fermact(), Inputparameters::tlsym);
	BOOST_REQUIRE_EQUAL(params.get_fermact_mp(), Inputparameters::iwasaki);
	BOOST_REQUIRE_EQUAL(params.get_kappa(), 0.25);
	BOOST_REQUIRE_EQUAL(params.get_mu(), 0.06);
	BOOST_REQUIRE_EQUAL(params.get_csw(), 0.1);
	BOOST_REQUIRE_EQUAL(params.get_iter0(), 2);
	BOOST_REQUIRE_EQUAL(params.get_iter1(), 3);
	BOOST_REQUIRE_EQUAL(params.get_kappa_mp(), 0.1125);
	BOOST_REQUIRE_EQUAL(params.get_mu_mp(), 0.016);
	BOOST_REQUIRE_EQUAL(params.get_csw_mp(), 0.2);
	BOOST_REQUIRE_EQUAL(params.get_iter0_mp(), 4);
	BOOST_REQUIRE_EQUAL(params.get_iter1_mp(), 5);
	BOOST_REQUIRE_EQUAL(params.get_cgmax(), 100);
	BOOST_REQUIRE_EQUAL(params.get_cgmax_mp(), 103);
	BOOST_REQUIRE_EQUAL(params.get_theta_fermion_spatial(), 0.8);
	BOOST_REQUIRE_EQUAL(params.get_theta_fermion_temporal(), 0.9);
	BOOST_REQUIRE_EQUAL(params.get_chem_pot_re(), 0.11);
	BOOST_REQUIRE_EQUAL(params.get_chem_pot_im(), 0.12);
	BOOST_REQUIRE_EQUAL(params.get_use_eo(), false);
	//at the moment, only 2 solvers are implemented..
	BOOST_REQUIRE_EQUAL(params.get_use_cg(), true);
	BOOST_REQUIRE_EQUAL(params.get_use_cg_mp(), false);
	BOOST_REQUIRE_EQUAL(params.get_use_bicgstab_save(), true);
	BOOST_REQUIRE_EQUAL(params.get_use_bicgstab_save_mp(), false);
	BOOST_REQUIRE_EQUAL(params.get_use_pointsource(), false);
	BOOST_REQUIRE_EQUAL(params.get_use_gauge_only(), true);
	BOOST_REQUIRE_EQUAL(params.get_num_sources(), 3);
	BOOST_REQUIRE_EQUAL(params.get_pointsource_x(), 1);
	BOOST_REQUIRE_EQUAL(params.get_pointsource_y(), 2);
	BOOST_REQUIRE_EQUAL(params.get_pointsource_z(), 3);
	BOOST_REQUIRE_EQUAL(params.get_pointsource_t(), 4);
	BOOST_REQUIRE_EQUAL(params.get_solver_prec(), 1e-21);
	BOOST_REQUIRE_EQUAL(params.get_force_prec(), 1e-16);
	BOOST_REQUIRE_EQUAL(params.get_iter_refresh(), 105);
	BOOST_REQUIRE_EQUAL(params.get_iter_refresh_mp(), 107);

	//HMC specific parameters
	BOOST_REQUIRE_EQUAL(params.get_tau(), 0.51234);
	BOOST_REQUIRE_EQUAL(params.get_reversibility_check(), true);
	BOOST_REQUIRE_EQUAL(params.get_integrationsteps(0), 103);
	BOOST_REQUIRE_EQUAL(params.get_integrationsteps(1), 5);
	BOOST_REQUIRE_EQUAL(params.get_integrationsteps(2), 10);
	BOOST_REQUIRE_EQUAL(params.get_hmcsteps(), 98);
	BOOST_REQUIRE_EQUAL(params.get_num_timescales(), 3);
	BOOST_REQUIRE_EQUAL(params.get_integrator(0), Inputparameters::twomn);
	BOOST_REQUIRE_EQUAL(params.get_integrator(1), Inputparameters::leapfrog);
	BOOST_REQUIRE_EQUAL(params.get_integrator(2), Inputparameters::twomn);
	//this is the optimal value...
	BOOST_REQUIRE_EQUAL(params.get_lambda(0), 0.11931833275037836);
	BOOST_REQUIRE_EQUAL(params.get_lambda(1), 0.21931833275037836);
	BOOST_REQUIRE_EQUAL(params.get_lambda(2), 0.31931833275037836);

	//direction for the correlator
	BOOST_REQUIRE_EQUAL(params.get_corr_dir(), 2);

	BOOST_REQUIRE_EQUAL(params.get_use_same_rnd_numbers(), true);
	BOOST_REQUIRE_EQUAL(params.get_profile_solver(), true);

	BOOST_REQUIRE_EQUAL(params.get_use_merge_kernels_spinor(), true);
	BOOST_REQUIRE_EQUAL(params.get_use_merge_kernels_fermion(), true);
}

BOOST_AUTO_TEST_CASE(command_line1)
{
	const char* _params[] = {"foo", "--use_aniso=1", "--use_mp=true", "--use_smearing=false", "--nspace=32", "--ntime=12"};
	Inputparameters params(6, _params);
	BOOST_REQUIRE_EQUAL(params.get_use_aniso(), true);
	BOOST_REQUIRE_EQUAL(params.get_use_mp(), true);
	BOOST_REQUIRE_EQUAL(params.get_use_smearing(), false);
	BOOST_REQUIRE_EQUAL(params.get_nspace(), 32);
	BOOST_REQUIRE_EQUAL(params.get_ntime(), 12);
}

BOOST_AUTO_TEST_CASE(command_line2)
{
	const char* _params[] = {"foo", "--integrator0=foo"};
	BOOST_REQUIRE_THROW(Inputparameters(2, _params), Inputparameters::parse_aborted);
}

BOOST_AUTO_TEST_CASE(command_line3)
{
	const char* _params[] = {"foo", "--startcondition=foo"};
	BOOST_REQUIRE_THROW(Inputparameters(2, _params), Inputparameters::parse_aborted);
}

BOOST_AUTO_TEST_CASE(command_line4)
{
	const char* _params[] = {"foo", "--fermact=foo"};
	BOOST_REQUIRE_THROW(Inputparameters(2, _params), Inputparameters::parse_aborted);
}

BOOST_AUTO_TEST_CASE(command_line5)
{
	const char* _params[] = {"foo", "--disable-ocl-compiler-opt=true"};
	Inputparameters params(2, _params);
	BOOST_REQUIRE_EQUAL(params.is_ocl_compiler_opt_disabled(), true);
}

BOOST_AUTO_TEST_CASE(help)
{
	const char* _params[] = {"foo", "--help"};
	BOOST_REQUIRE_THROW(Inputparameters(2, _params), Inputparameters::parse_aborted);
}
