/** @file
 * Testcases for the Inputparameters class
 *
 * Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
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

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE meta::Inputparameters
#include <boost/test/unit_test.hpp>

#include "inputparameters.hpp"

#include "../common_header_files/types.h"

using namespace meta;

void checkDefaults(const Inputparameters & params)
{
	BOOST_REQUIRE_EQUAL(params.get_precision(), sizeof(double) * 8);

	BOOST_REQUIRE_EQUAL(params.get_selected_devices().size(), 0);
	BOOST_REQUIRE_EQUAL(params.get_device_count(), 0);
	BOOST_REQUIRE_EQUAL(params.get_use_gpu(), true);
	BOOST_REQUIRE_EQUAL(params.get_use_cpu(), true);
	BOOST_REQUIRE_EQUAL(params.get_enable_profiling(), false);

	BOOST_REQUIRE_EQUAL(params.get_use_aniso(), false);
	BOOST_REQUIRE_EQUAL(params.get_use_chem_pot_re(), false);
	BOOST_REQUIRE_EQUAL(params.get_use_chem_pot_im(), false);
	BOOST_REQUIRE_EQUAL(params.get_use_smearing(), false);
	BOOST_REQUIRE_EQUAL(params.get_use_mp(), false);
	BOOST_REQUIRE_EQUAL(params.get_nspace(), 4);
	BOOST_REQUIRE_EQUAL(params.get_ntime(), 8);

	BOOST_REQUIRE_EQUAL(params.get_startcondition(), common::cold_start);
	BOOST_REQUIRE_EQUAL(params.get_writefrequency(), 1);
	BOOST_REQUIRE_EQUAL(params.get_savefrequency(), 100);
	BOOST_REQUIRE_EQUAL(params.get_sourcefile(), "conf.00000");
	BOOST_REQUIRE_EQUAL(params.get_ignore_checksum_errors(), false);
	BOOST_REQUIRE_EQUAL(params.get_print_to_screen(), false);
	//This is obvious!!!
	BOOST_REQUIRE_EQUAL(params.get_host_seed(), 4815);

	//gaugefield parameters
	BOOST_REQUIRE_EQUAL(params.get_beta(), 4.0);
	BOOST_REQUIRE_EQUAL(params.get_rho(), 0.);
	BOOST_REQUIRE_EQUAL(params.get_rho_iter(), 0);
	BOOST_REQUIRE_EQUAL(params.get_gaugeact(), common::action::wilson);

	//heatbath parameters
	BOOST_REQUIRE_EQUAL(params.get_thermalizationsteps(), 0);
	BOOST_REQUIRE_EQUAL(params.get_heatbathsteps(), 1000);
	BOOST_REQUIRE_EQUAL(params.get_overrelaxsteps(), 1);
	BOOST_REQUIRE_EQUAL(params.get_xi(), 1);
	BOOST_REQUIRE_EQUAL(params.get_measure_transportcoefficient_kappa(), false);
	BOOST_REQUIRE_EQUAL(params.get_measure_rectangles(), false);

	//fermionic parameters
	BOOST_REQUIRE_EQUAL(params.get_fermact(), common::action::wilson);
	BOOST_REQUIRE_EQUAL(params.get_fermact_mp(), common::action::wilson);
	BOOST_REQUIRE_EQUAL(params.get_kappa(), 0.125);
	BOOST_REQUIRE_EQUAL(params.get_mu(), 0.006);
	BOOST_REQUIRE_EQUAL(params.get_csw(), 0.);
	BOOST_REQUIRE_EQUAL(params.get_kappa_mp(), 0.125);
	BOOST_REQUIRE_EQUAL(params.get_mu_mp(), 0.006);
	BOOST_REQUIRE_EQUAL(params.get_csw_mp(), 0.);
	BOOST_REQUIRE_EQUAL(params.get_cgmax(), 1000);
	BOOST_REQUIRE_EQUAL(params.get_cgmax_mp(), 1000);
	BOOST_REQUIRE_EQUAL(params.get_theta_fermion_spatial(), 0.);
	BOOST_REQUIRE_EQUAL(params.get_theta_fermion_temporal(), 0.);
	BOOST_REQUIRE_EQUAL(params.get_chem_pot_re(), 0.);
	BOOST_REQUIRE_EQUAL(params.get_chem_pot_im(), 0.);
	BOOST_REQUIRE_EQUAL(params.get_use_eo(), true);
	//at the moment, only 2 solvers are implemented..
	BOOST_REQUIRE_EQUAL(params.get_solver(), common::bicgstab);
	BOOST_REQUIRE_EQUAL(params.get_use_gauge_only(), false);
	BOOST_REQUIRE_EQUAL(params.get_num_sources(), 12);
	BOOST_REQUIRE_EQUAL(params.get_source_x(), 0);
	BOOST_REQUIRE_EQUAL(params.get_source_y(), 0);
	BOOST_REQUIRE_EQUAL(params.get_source_z(), 0);
	BOOST_REQUIRE_EQUAL(params.get_source_t(), 0);
#ifdef _USEDOUBLEPREC_
	BOOST_REQUIRE_EQUAL(params.get_solver_prec(), 1e-23);
	BOOST_REQUIRE_EQUAL(params.get_force_prec(), 1e-12);
#else
	BOOST_REQUIRE_EQUAL(params.get_solver_prec(), 1e-16);
	BOOST_REQUIRE_EQUAL(params.get_force_prec(), 1e-8);
#endif
	BOOST_REQUIRE_EQUAL(params.get_iter_refresh(), 100);
	BOOST_REQUIRE_EQUAL(params.get_iter_refresh_mp(), 100);
	BOOST_REQUIRE_EQUAL(params.get_benchmarksteps(), 500);

	//HMC specific parameters
	BOOST_REQUIRE_EQUAL(params.get_tau(), 0.5);
	BOOST_REQUIRE_EQUAL(params.get_reversibility_check(), false);
	BOOST_REQUIRE_EQUAL(params.get_integrationsteps(0), 10);
	BOOST_REQUIRE_EQUAL(params.get_integrationsteps(1), 10);
	BOOST_REQUIRE_EQUAL(params.get_integrationsteps(2), 10);
	BOOST_REQUIRE_EQUAL(params.get_hmcsteps(), 10);
	BOOST_REQUIRE_EQUAL(params.get_num_timescales(), 1);
	BOOST_REQUIRE_EQUAL(params.get_integrator(0), common::leapfrog);
	BOOST_REQUIRE_EQUAL(params.get_integrator(1), common::leapfrog);
	BOOST_REQUIRE_EQUAL(params.get_integrator(2), common::leapfrog);
	BOOST_REQUIRE_EQUAL(params.get_lambda(0), 0.1931833275037836); //this is the optimal value...
	BOOST_REQUIRE_EQUAL(params.get_lambda(1), 0.1931833275037836);
	BOOST_REQUIRE_EQUAL(params.get_lambda(2), 0.1931833275037836);

	//RHMC specific parameters
	BOOST_REQUIRE_EQUAL(params.get_md_approx_ord(), 8);
	BOOST_REQUIRE_EQUAL(params.get_metro_approx_ord(), 15);
	BOOST_REQUIRE_EQUAL(params.get_findminmax_iteration_block_size(), 25);
	BOOST_REQUIRE_EQUAL(params.get_findminmax_max(), 5000);
	BOOST_REQUIRE_EQUAL(params.get_findminmax_prec(), 0.001);
	BOOST_REQUIRE_EQUAL(params.get_conservative(), false);
	BOOST_REQUIRE_EQUAL(params.get_num_tastes(), 2);
	BOOST_REQUIRE_EQUAL(params.get_num_tastes_decimal_digits(), 0);
	BOOST_REQUIRE_EQUAL(params.get_num_pseudofermions(), 1);
	BOOST_REQUIRE_EQUAL(params.get_approx_lower(), 1.e-5);
	BOOST_REQUIRE_EQUAL(params.get_approx_upper(), 1.);
	BOOST_REQUIRE_EQUAL(params.get_rhmcsteps(), 10);
	BOOST_REQUIRE_EQUAL(params.get_approx_heatbath_file(), "Approx_Heatbath");
	BOOST_REQUIRE_EQUAL(params.get_approx_md_file(), "Approx_MD");
	BOOST_REQUIRE_EQUAL(params.get_approx_metropolis_file(), "Approx_Metropolis");
	BOOST_REQUIRE_EQUAL(params.get_read_rational_approximations_from_file(), true);

	//direction for the correlator
	BOOST_REQUIRE_EQUAL(params.get_corr_dir(), 3);

	BOOST_REQUIRE_EQUAL(params.get_use_same_rnd_numbers(), false);
	BOOST_REQUIRE_EQUAL(params.get_profile_solver(), false);
	BOOST_REQUIRE_EQUAL(params.get_test_ref_value(), 0.);
	BOOST_REQUIRE_EQUAL(params.get_test_ref_value2(), 0.);

	BOOST_REQUIRE_EQUAL(params.is_ocl_compiler_opt_disabled(), false);

	BOOST_REQUIRE_EQUAL(params.get_use_merge_kernels_spinor(), false);
	BOOST_REQUIRE_EQUAL(params.get_use_merge_kernels_fermion(), false);
	BOOST_REQUIRE_EQUAL(params.get_use_rec12(), false);

	BOOST_REQUIRE_EQUAL(params.get_log_level(), "ALL");
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
	BOOST_REQUIRE_EQUAL(params.get_use_gpu(), false);
	BOOST_REQUIRE_EQUAL(params.get_use_cpu(), false);
	BOOST_REQUIRE_EQUAL(params.get_enable_profiling(), true);

	BOOST_REQUIRE_EQUAL(params.get_use_aniso(), true);
	BOOST_REQUIRE_EQUAL(params.get_use_chem_pot_re(), true);
	BOOST_REQUIRE_EQUAL(params.get_use_chem_pot_im(), true);
	BOOST_REQUIRE_EQUAL(params.get_use_smearing(), true);
	BOOST_REQUIRE_EQUAL(params.get_use_mp(), true);
	BOOST_REQUIRE_EQUAL(params.get_nspace(), 32);
	BOOST_REQUIRE_EQUAL(params.get_ntime(), 12);

	BOOST_REQUIRE_EQUAL(params.get_startcondition(), common::start_from_source);
	BOOST_REQUIRE_EQUAL(params.get_writefrequency(), 10);
	BOOST_REQUIRE_EQUAL(params.get_savefrequency(), 10);
	BOOST_REQUIRE_EQUAL(params.get_sourcefile(), "conf.11111");
	BOOST_REQUIRE_EQUAL(params.get_ignore_checksum_errors(), true);
	BOOST_REQUIRE_EQUAL(params.get_print_to_screen(), true);
	//This is obvious!!!
	BOOST_REQUIRE_EQUAL(params.get_host_seed(), 42);

	//gaugefield parameters
	BOOST_REQUIRE_EQUAL(params.get_beta(), 4.1);
	BOOST_REQUIRE_EQUAL(params.get_rho(), 0.1);
	BOOST_REQUIRE_EQUAL(params.get_rho_iter(), 4);
	BOOST_REQUIRE_EQUAL(params.get_gaugeact(), common::action::twistedmass);

	//heatbath parameters
	BOOST_REQUIRE_EQUAL(params.get_thermalizationsteps(), 10);
	BOOST_REQUIRE_EQUAL(params.get_heatbathsteps(), 100);
	BOOST_REQUIRE_EQUAL(params.get_overrelaxsteps(), 10);
	BOOST_REQUIRE_EQUAL(params.get_xi(), 2);
	BOOST_REQUIRE_EQUAL(params.get_measure_transportcoefficient_kappa(), true);

	//fermionic parameters
	BOOST_REQUIRE_EQUAL(params.get_fermact(), common::action::tlsym);
	BOOST_REQUIRE_EQUAL(params.get_fermact_mp(), common::action::iwasaki);
	BOOST_REQUIRE_EQUAL(params.get_kappa(), 0.25);
	BOOST_REQUIRE_EQUAL(params.get_mu(), 0.06);
	BOOST_REQUIRE_EQUAL(params.get_csw(), 0.1);
	BOOST_REQUIRE_EQUAL(params.get_kappa_mp(), 0.1125);
	BOOST_REQUIRE_EQUAL(params.get_mu_mp(), 0.016);
	BOOST_REQUIRE_EQUAL(params.get_csw_mp(), 0.2);
	BOOST_REQUIRE_EQUAL(params.get_cgmax(), 100);
	BOOST_REQUIRE_EQUAL(params.get_cgmax_mp(), 103);
	BOOST_REQUIRE_EQUAL(params.get_theta_fermion_spatial(), 0.8);
	BOOST_REQUIRE_EQUAL(params.get_theta_fermion_temporal(), 0.9);
	BOOST_REQUIRE_EQUAL(params.get_chem_pot_re(), 0.11);
	BOOST_REQUIRE_EQUAL(params.get_chem_pot_im(), 0.12);
	BOOST_REQUIRE_EQUAL(params.get_use_eo(), false);
	//at the moment, only 2 solvers are implemented..
	BOOST_REQUIRE_EQUAL(params.get_solver(), common::cg);
	BOOST_REQUIRE_EQUAL(params.get_use_gauge_only(), true);
	BOOST_REQUIRE_EQUAL(params.get_num_sources(), 3);
	BOOST_REQUIRE_EQUAL(params.get_source_x(), 1);
	BOOST_REQUIRE_EQUAL(params.get_source_y(), 2);
	BOOST_REQUIRE_EQUAL(params.get_source_z(), 3);
	BOOST_REQUIRE_EQUAL(params.get_source_t(), 4);
	BOOST_REQUIRE_EQUAL(params.get_solver_prec(), 1e-21);
	BOOST_REQUIRE_EQUAL(params.get_force_prec(), 1e-16);
	BOOST_REQUIRE_EQUAL(params.get_iter_refresh(), 105);
	BOOST_REQUIRE_EQUAL(params.get_iter_refresh_mp(), 107);

	BOOST_REQUIRE_EQUAL(params.get_benchmarksteps(), 10);

	//HMC specific parameters
	BOOST_REQUIRE_EQUAL(params.get_tau(), 0.51234);
	BOOST_REQUIRE_EQUAL(params.get_reversibility_check(), true);
	BOOST_REQUIRE_EQUAL(params.get_integrationsteps(0), 103);
	BOOST_REQUIRE_EQUAL(params.get_integrationsteps(1), 5);
	BOOST_REQUIRE_EQUAL(params.get_integrationsteps(2), 10);
	BOOST_REQUIRE_EQUAL(params.get_hmcsteps(), 98);
	BOOST_REQUIRE_EQUAL(params.get_num_timescales(), 3);
	BOOST_REQUIRE_EQUAL(params.get_integrator(0), common::twomn);
	BOOST_REQUIRE_EQUAL(params.get_integrator(1), common::leapfrog);
	BOOST_REQUIRE_EQUAL(params.get_integrator(2), common::twomn);
	//this is the optimal value...
	BOOST_REQUIRE_EQUAL(params.get_lambda(0), 0.11931833275037836);
	BOOST_REQUIRE_EQUAL(params.get_lambda(1), 0.21931833275037836);
	BOOST_REQUIRE_EQUAL(params.get_lambda(2), 0.31931833275037836);

	//direction for the correlator
	BOOST_REQUIRE_EQUAL(params.get_corr_dir(), 2);

	BOOST_REQUIRE_EQUAL(params.get_use_same_rnd_numbers(), true);
	BOOST_REQUIRE_EQUAL(params.get_profile_solver(), true);
	BOOST_REQUIRE_EQUAL(params.get_test_ref_value(), 1.337);
	BOOST_REQUIRE_EQUAL(params.get_test_ref_value2(), 2.448);

	BOOST_REQUIRE_EQUAL(params.get_use_merge_kernels_spinor(), true);
	BOOST_REQUIRE_EQUAL(params.get_use_merge_kernels_fermion(), true);
	BOOST_REQUIRE_EQUAL(params.get_use_rec12(), true);
}

BOOST_AUTO_TEST_CASE(command_line1)
{
	const char* _params[] = {"foo", "--use_aniso=1", "--use_mp=true", "--use_smearing=false", "--nspace=32", "--ntime=12", "--use_gpu=false", "--use_cpu=false"};
	Inputparameters params(8, _params);
	BOOST_REQUIRE_EQUAL(params.get_use_aniso(), true);
	BOOST_REQUIRE_EQUAL(params.get_use_mp(), true);
	BOOST_REQUIRE_EQUAL(params.get_use_smearing(), false);
	BOOST_REQUIRE_EQUAL(params.get_nspace(), 32);
	BOOST_REQUIRE_EQUAL(params.get_ntime(), 12);
	BOOST_REQUIRE_EQUAL(params.get_use_gpu(), false);
	BOOST_REQUIRE_EQUAL(params.get_use_cpu(), false);
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

BOOST_AUTO_TEST_CASE(command_line6)
{
	const char* _params[] = {"foo", "--log-level=ERROR"};
	Inputparameters params(2, _params);
	BOOST_REQUIRE_EQUAL(params.get_log_level(), "ERROR");
}

BOOST_AUTO_TEST_CASE(aliases)
{
	const char* _params[] = {"foo", "test_input_aliases"};
	Inputparameters params(2, _params);
	BOOST_REQUIRE_EQUAL(params.get_nspace(), 12);
	BOOST_REQUIRE_EQUAL(params.get_ntime(), 32);
	BOOST_REQUIRE_EQUAL(params.get_use_eo(), false);
	BOOST_REQUIRE_EQUAL(params.get_thermalizationsteps(), 13);
	BOOST_REQUIRE_EQUAL(params.get_use_gauge_only(), true);
	BOOST_REQUIRE_EQUAL(params.get_test_ref_value(), 42.);
	BOOST_REQUIRE_EQUAL(params.get_test_ref_value2(), 53.);
	BOOST_REQUIRE_EQUAL(params.get_theta_fermion_spatial(), 3.);
	BOOST_REQUIRE_EQUAL(params.get_theta_fermion_temporal(), 1.5);
}

BOOST_AUTO_TEST_CASE(help)
{
	const char* _params[] = {"foo", "--help"};
	BOOST_REQUIRE_THROW(Inputparameters(2, _params), Inputparameters::parse_aborted);
}

BOOST_AUTO_TEST_CASE(hmcParameters)
{
	const char* _params[] = {"foo", "--help"};
	BOOST_REQUIRE_THROW(Inputparameters(2, _params, "hmc"), Inputparameters::parse_aborted);
}

BOOST_AUTO_TEST_CASE(rhmcParameters)
{
	const char* _params[] = {"foo", "--help"};
	BOOST_REQUIRE_THROW(Inputparameters(2, _params, "rhmc"), Inputparameters::parse_aborted);
}

BOOST_AUTO_TEST_CASE(heatbathParameters)
{
	const char* _params[] = {"foo", "--help"};
	BOOST_REQUIRE_THROW(Inputparameters(2, _params, "su3heatbath"), Inputparameters::parse_aborted);
}

BOOST_AUTO_TEST_CASE(inverterParameters)
{
	const char* _params[] = {"foo", "--help"};
	BOOST_REQUIRE_THROW(Inputparameters(2, _params, "inverter"), Inputparameters::parse_aborted);
}

BOOST_AUTO_TEST_CASE(gaugeobservablesParameters)
{
	const char* _params[] = {"foo", "--help"};
	BOOST_REQUIRE_THROW(Inputparameters(2, _params, "gaugeobservables"), Inputparameters::parse_aborted);
}
