/** @file
 * Testcases for the Inputparameters class
 *
 * Copyright (c) 2014 Christopher Pinke
 * Copyright (c) 2014 Matthias Bach
 * Copyright (c) 2018 Alessandro Sciarra
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with CL2QCD. If not, see <http://www.gnu.org/licenses/>.
 */

// use the boost test framework
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE meta::parametersHmc
#include <boost/test/unit_test.hpp>

#include "../common_header_files/types.hpp"

#include "parametersHmc.hpp"

using namespace meta;

void checkDefaults(const ParametersHmc params)
{
	BOOST_REQUIRE_EQUAL(params.get_tau(), 0.5);
	BOOST_REQUIRE_EQUAL(params.get_reversibility_check(), false);
	BOOST_REQUIRE_EQUAL(params.get_integrationsteps(0), 10);
	BOOST_REQUIRE_EQUAL(params.get_integrationsteps(1), 10);
	BOOST_REQUIRE_EQUAL(params.get_integrationsteps(2), 10);
	BOOST_REQUIRE_EQUAL(params.get_hmcsteps(), 10);
	BOOST_REQUIRE_EQUAL(params.get_num_timescales(), 1);
	BOOST_REQUIRE_EQUAL(params.get_integrator(0), ParametersHmc::leapfrog);
	BOOST_REQUIRE_EQUAL(params.get_integrator(1), ParametersHmc::leapfrog);
	BOOST_REQUIRE_EQUAL(params.get_integrator(2), ParametersHmc::leapfrog);
	//this is the optimal value...
	BOOST_REQUIRE_EQUAL(params.get_lambda(0), 0.1931833275037836);
	BOOST_REQUIRE_EQUAL(params.get_lambda(1), 0.1931833275037836);
	BOOST_REQUIRE_EQUAL(params.get_lambda(2), 0.1931833275037836);
}

BOOST_AUTO_TEST_CASE(defaults)
{
	const char* _params[] = {"foo"};
	ParametersHmc params(1, _params);
	checkDefaults(params);
}

BOOST_AUTO_TEST_CASE(input_file1)
{
	// empty input file
	const char* _params[] = {"foo", "test_input_1"};
	ParametersHmc params(2, _params);
	checkDefaults(params);
}

// clang-format off
/*
BOOST_AUTO_TEST_CASE(input_file2)
{
  // filled input file
  const char* _params2[] = {"foo", "test_input_2"};
  ParametersHmc params(2, _params2);
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

  BOOST_REQUIRE_EQUAL(params.get_startcondition(), Inputparameters::start_from_source);
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
  BOOST_REQUIRE_EQUAL(params.get_gaugeact(), Inputparameters::twistedmass);

  //heatbath parameters
  BOOST_REQUIRE_EQUAL(params.get_thermalizationsteps(), 10);
  BOOST_REQUIRE_EQUAL(params.get_heatbathsteps(), 100);
  BOOST_REQUIRE_EQUAL(params.get_overrelaxsteps(), 10);
  BOOST_REQUIRE_EQUAL(params.get_xi(), 2);
  BOOST_REQUIRE_EQUAL(params.get_measure_transportcoefficient_kappa(), true);

  //fermionic parameters
  BOOST_REQUIRE_EQUAL(params.get_fermact(), Inputparameters::tlsym);
  BOOST_REQUIRE_EQUAL(params.get_fermact_mp(), Inputparameters::iwasaki);
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
  BOOST_REQUIRE_EQUAL(params.get_solver(), Inputparameters::cg);
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
  BOOST_REQUIRE_EQUAL(params.get_test_ref_value(), 1.337);
  BOOST_REQUIRE_EQUAL(params.get_test_ref_value2(), 2.448);

  BOOST_REQUIRE_EQUAL(params.get_use_merge_kernels_spinor(), true);
  BOOST_REQUIRE_EQUAL(params.get_use_merge_kernels_fermion(), true);
  BOOST_REQUIRE_EQUAL(params.get_use_rec12(), true);
}

BOOST_AUTO_TEST_CASE(command_line1)
{
  const char* _params[] = {"foo", "--use_aniso=1", "--use_mp=true", "--use_smearing=false", "--nspace=32", "--ntime=12", "--use_gpu=false", "--use_cpu=false"};
  ParametersHmc params(8, _params);
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
*/
// clang-format on
