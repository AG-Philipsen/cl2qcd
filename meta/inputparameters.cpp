/** @file
 * Input file handling implementation
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

#include "inputparameters.hpp"

#include <map>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
namespace po = boost::program_options;

#include "config_file_normalizer.hpp"

using namespace meta;

/**
 * Get the action describe by the given string.
 */
static Inputparameters::action get_action(std::string);
/**
 * Get the integrator describe by the given string.
 */
static Inputparameters::integrator get_integrator(std::string);
/**
 * Get the startcondition describe by the given string.
 */
static Inputparameters::startcondition get_startcondition(std::string);
/**
 * Get the solver describe by the given string.
 */
static Inputparameters::solver get_solver(std::string);
/**
 * Get the sourcetype given in string
 */
static Inputparameters::sourcetypes get_sourcetype(std::string s);
/**
 * Get the pbp_version given in string.
 */
static Inputparameters::pbp_version get_pbp_version(std::string s);
/**
 * Get the sourcecontent given in string.
 */
static Inputparameters::sourcecontents get_sourcecontent(std::string s);
/**
 * Adds all alternative option names to the ConfigFileNormlizer instance
 */
static void add_option_aliases(meta::ConfigFileNormalizer * const);

size_t Inputparameters::get_precision() const noexcept
{
  return precision;
}

const std::vector<int> Inputparameters::get_selected_devices() const noexcept
{
  return selected_devices;
}
int Inputparameters::get_device_count() const noexcept
{
  return device_count;
}
bool Inputparameters::get_use_gpu() const noexcept
{
  return use_gpu;
}
bool Inputparameters::get_use_cpu() const noexcept
{
  return use_cpu;
}

bool Inputparameters::get_use_aniso() const noexcept
{
  return use_aniso;
}
bool Inputparameters::get_use_chem_pot_re() const noexcept
{
  return use_chem_pot_re;
}
bool Inputparameters::get_use_chem_pot_im() const noexcept
{
  return use_chem_pot_im;
}
bool Inputparameters::get_use_smearing() const noexcept
{
  return use_smearing;
}
bool Inputparameters::get_use_mp() const noexcept
{
  return use_mp;
}
int Inputparameters::get_nspace() const noexcept
{
  return nspace;
}
int Inputparameters::get_ntime() const noexcept
{
  return ntime;
}

Inputparameters::startcondition Inputparameters::get_startcondition() const noexcept
{
  return _startcondition;
}
int Inputparameters::get_writefrequency() const noexcept
{
  return writefrequency;
}
int Inputparameters::get_savefrequency() const noexcept
{
  return savefrequency;
}
std::string Inputparameters::get_sourcefile() const noexcept
{
  return sourcefile;
}
bool Inputparameters::get_ignore_checksum_errors() const noexcept
{
	return ignore_checksum_errors;
}
bool Inputparameters::get_print_to_screen() const noexcept
{
  return print_to_screen;
}
uint32_t Inputparameters::get_host_seed() const noexcept
{
  return host_seed;
}
std::string Inputparameters::get_initial_prng_state() const noexcept
{
  return initial_prng_state;
}

//gaugefield parameters
double Inputparameters::get_beta() const noexcept
{
  return beta;
}
double Inputparameters::get_rho() const noexcept
{
  return rho;
}
int Inputparameters::get_rho_iter() const noexcept
{
  return rho_iter;
}
Inputparameters::action Inputparameters::get_gaugeact() const noexcept
{
  return gaugeact;
}

//heatbath parameters
int Inputparameters::get_thermalizationsteps() const noexcept
{
  return thermalizationsteps;
}
int Inputparameters::get_heatbathsteps() const noexcept
{
  return heatbathsteps;
}
int Inputparameters::get_overrelaxsteps() const noexcept
{
  return overrelaxsteps;
}
int Inputparameters::get_xi() const noexcept
{
  return xi;
}
bool Inputparameters::get_measure_transportcoefficient_kappa() const noexcept
{
	return measure_transportcoefficient_kappa;
}

//fermionic parameters
Inputparameters::action Inputparameters::get_fermact() const noexcept
{
  return fermact;
}
Inputparameters::action Inputparameters::get_fermact_mp() const noexcept
{
  return fermact_mp;
}
double Inputparameters::get_kappa() const noexcept
{
  return kappa;
}
double Inputparameters::get_mass() const noexcept
{
  return mass;
}
double Inputparameters::get_mu() const noexcept
{
  return mu;
}
double Inputparameters::get_csw() const noexcept
{
  return csw;
}
double Inputparameters::get_kappa_mp() const noexcept
{
  return kappa_mp;
}
double Inputparameters::get_mu_mp() const noexcept
{
  return mu_mp;
}
double Inputparameters::get_csw_mp() const noexcept
{
  return csw_mp;
}
int Inputparameters::get_cgmax() const noexcept
{
  return cgmax;
}
int Inputparameters::get_cgmax_mp() const noexcept
{
  return cgmax_mp;
}
double Inputparameters::get_theta_fermion_spatial() const noexcept
{
  return theta_fermion_spatial;
}
double Inputparameters::get_theta_fermion_temporal() const noexcept
{
  return theta_fermion_temporal;
}
double Inputparameters::get_chem_pot_re() const noexcept
{
  return chem_pot_re;
}
double Inputparameters::get_chem_pot_im() const noexcept
{
  return chem_pot_im;
}
bool Inputparameters::get_use_eo() const noexcept
{
  return use_eo;
}
bool Inputparameters::get_use_gauge_only() const noexcept
{
  return use_gauge_only;
}
int Inputparameters::get_num_sources() const noexcept
{
  return num_sources;
}
int Inputparameters::get_source_x() const noexcept
{
  return source_x;
}
int Inputparameters::get_source_y() const noexcept
{
  return source_y;
}
int Inputparameters::get_source_z() const noexcept
{
  return source_z;
}
int Inputparameters::get_source_t() const noexcept
{
  return source_t;
}
bool Inputparameters::get_place_sources_on_host() const noexcept
{
	return place_sources_on_host;
}

double Inputparameters::get_solver_prec() const noexcept
{
  return solver_prec;
}
double Inputparameters::get_force_prec() const noexcept
{
  return force_prec;
}
int Inputparameters::get_iter_refresh() const noexcept
{
  return iter_refresh;
}
int Inputparameters::get_iter_refresh_mp() const noexcept
{
  return iter_refresh_mp;
}
Inputparameters::solver Inputparameters::get_solver() const noexcept
{
  return _solver;
}
Inputparameters::solver Inputparameters::get_solver_mp() const noexcept
{
  return _solver_mp;
}

//HMC specific parameters
double Inputparameters::get_tau() const noexcept
{
  return tau;
}
bool Inputparameters::get_reversibility_check() const noexcept
{
  return reversibility_check;
}
int Inputparameters::get_integrationsteps(size_t timescale) const noexcept
{
switch(timescale) {
case 0:
return integrationsteps0;
case 1:
return integrationsteps1;
case 2:
return integrationsteps2;
default:
throw std::out_of_range("No such timescale");
}
}
int Inputparameters::get_hmcsteps() const noexcept
{
  return hmcsteps;
}
int Inputparameters::get_num_timescales() const noexcept
{
  return num_timescales;
}
Inputparameters::integrator Inputparameters::get_integrator(size_t timescale) const noexcept
{
switch(timescale) {
case 0:
return integrator0;
case 1:
return integrator1;
case 2:
return integrator2;
default:
throw std::out_of_range("No such timescale");
}
}
double Inputparameters::get_lambda(size_t timescale) const noexcept
{
switch(timescale) {
case 0:
return lambda0;
case 1:
return lambda1;
case 2:
return lambda2;
default:
throw std::out_of_range("No such timescale");
}
}

//RHMC specific parameters
int Inputparameters::get_md_approx_ord() const noexcept
{
  return md_approx_ord;
}
int Inputparameters::get_metro_approx_ord() const noexcept
{
  return metro_approx_ord;
}
int Inputparameters::get_findminmax_iteration_block_size() const noexcept
{
  return findminmax_iteration_block_size;
}
int Inputparameters::get_findminmax_max() const noexcept
{
  return findminmax_max;
}
double Inputparameters::get_findminmax_prec() const noexcept
{
  return findminmax_prec;
}
bool Inputparameters::get_conservative() const noexcept
{
  return conservative;
}
int Inputparameters::get_num_tastes() const noexcept
{
  return num_tastes;
}
double Inputparameters::get_approx_lower() const noexcept
{
  return approx_lower;
}
double Inputparameters::get_approx_upper() const noexcept
{
  return approx_upper;
}
int Inputparameters::get_rhmcsteps() const noexcept
{
  return rhmcsteps;
}

//direction for the correlator
int Inputparameters::get_corr_dir() const noexcept
{
  return corr_dir;
}

bool Inputparameters::get_use_same_rnd_numbers() const noexcept
{
  return use_same_rnd_numbers;
}
bool Inputparameters::get_profile_solver() const noexcept
{
  return profile_solver;
}

bool Inputparameters::is_ocl_compiler_opt_disabled() const noexcept
{
  return ocl_compiler_opt_disabled;
}
double Inputparameters::get_test_ref_value() const noexcept
{
  return test_ref_value;
}

std::string Inputparameters::get_log_level() const noexcept
{
  return log_level;
}
bool Inputparameters::get_use_merge_kernels_fermion() const noexcept
{
  return use_merge_kernels_fermion;
}
bool Inputparameters::get_use_merge_kernels_spinor() const noexcept
{
  return use_merge_kernels_spinor;
}
bool Inputparameters::get_use_rec12() const noexcept
{
  return use_rec12;
}

//parameters to read in gauge configurations
bool Inputparameters::get_read_multiple_configs() const noexcept
{
  return read_multiple_configs;
}
int Inputparameters::get_config_read_start() const noexcept
{
  return config_read_start;
}
int Inputparameters::get_config_read_end() const noexcept
{
  return config_read_end;
}
int Inputparameters::get_config_read_incr() const noexcept
{
  return config_read_incr;
}
int Inputparameters::get_config_number_digits() const noexcept
{
  return config_number_digits;
}
std::string Inputparameters::get_prng_prefix() const noexcept
{
  return prng_prefix;
}
std::string Inputparameters::get_prng_postfix() const noexcept
{
  return prng_postfix;
}
std::string Inputparameters::get_config_prefix() const noexcept
{
  return config_prefix;
}
std::string Inputparameters::get_config_postfix() const noexcept
{
  return config_postfix;
}
std::string Inputparameters::get_ferm_obs_corr_prefix() const noexcept
{
  return ferm_obs_corr_prefix;
}
std::string Inputparameters::get_ferm_obs_corr_postfix() const noexcept
{
  return ferm_obs_corr_postfix;
}
std::string Inputparameters::get_ferm_obs_pbp_prefix() const noexcept
{
  return ferm_obs_pbp_prefix;
}
std::string Inputparameters::get_ferm_obs_pbp_postfix() const noexcept
{
  return ferm_obs_pbp_postfix;
}
std::string Inputparameters::get_gauge_obs_prefix() const noexcept
{
  return gauge_obs_prefix;
}
std::string Inputparameters::get_gauge_obs_postfix() const noexcept
{
  return gauge_obs_postfix;
}
bool Inputparameters::get_ferm_obs_to_single_file() const noexcept
{
  return ferm_obs_to_single_file;
}
bool Inputparameters::get_gauge_obs_to_single_file() const noexcept
{
  return gauge_obs_to_single_file;
}
std::string Inputparameters::get_hmc_obs_prefix() const noexcept
{
  return hmc_obs_prefix;
}
std::string Inputparameters::get_hmc_obs_postfix() const noexcept
{
  return hmc_obs_postfix;
}
bool Inputparameters::get_hmc_obs_to_single_file() const noexcept
{
  return hmc_obs_to_single_file;
}
std::string Inputparameters::get_rhmc_obs_prefix() const noexcept
{
  return rhmc_obs_prefix;
}
std::string Inputparameters::get_rhmc_obs_postfix() const noexcept
{
  return rhmc_obs_postfix;
}
bool Inputparameters::get_rhmc_obs_to_single_file() const noexcept
{
  return rhmc_obs_to_single_file;
}
Inputparameters::sourcetypes Inputparameters::get_sourcetype() const noexcept
{
  return sourcetype;
}
Inputparameters::sourcecontents Inputparameters::get_sourcecontent() const noexcept
{
  return sourcecontent;
}
bool Inputparameters::get_measure_correlators() const noexcept
{
  return measure_correlators;
}
bool Inputparameters::get_measure_pbp() const noexcept
{
  return measure_pbp;
}
Inputparameters::pbp_version Inputparameters::get_pbp_version() const noexcept
{
  return pbp_version_;
}
int Inputparameters::get_cg_iteration_block_size() const noexcept
{
  return cg_iteration_block_size;
}
bool Inputparameters::get_cg_use_async_copy() const noexcept
{
  return cg_use_async_copy;
}
int Inputparameters::get_cg_minimum_iteration_count() const noexcept
{
  return cg_minimum_iteration_count;
}

Inputparameters::Inputparameters(int argc, const char** argv)
{
	logger.info() << "read in parameters...";
	/**
	 * First handle all the stuff that can only be done on the cmd-line.
	 * We need that to get the option file.
	 */
	po::options_description cmd_opts("Generic options");
	cmd_opts.add_options()
	("help,h", "Produce this help message")
	("input-file", po::value<std::string>(), "File containing the input parameters");
	// TODO add log-level etc
	po::positional_options_description pos_opts;
	pos_opts.add("input-file", 1);

	po::options_description config("Configuration options");
	config.add_options()
	("prec", po::value<size_t>(&precision)->default_value(sizeof(double) * 8))

	("device,d", po::value<std::vector<int>>(&selected_devices), "ID of a divice to use. Can be specified multiple times.")
	("num_dev", po::value<int>(&device_count)->default_value(0), "Maximum number of devices to use.")
	("use_gpu", po::value<bool>(&use_gpu)->default_value(true), "Use GPUs")
	("use_cpu", po::value<bool>(&use_cpu)->default_value(true), "Use CPUs")

	("use_aniso", po::value<bool>(&use_aniso)->default_value(false))
	("use_chem_pot_re", po::value<bool>(&use_chem_pot_re)->default_value(false))
	("use_chem_pot_im", po::value<bool>(&use_chem_pot_im)->default_value(false))
	("use_smearing", po::value<bool>(&use_smearing)->default_value(false))
	("use_mp", po::value<bool>(&use_mp)->default_value(false))
	("nspace", po::value<int>(&nspace)->default_value(4))
	("ntime", po::value<int>(&ntime)->default_value(8))

	("startcondition", po::value<std::string>()->default_value("cold_start"))
	("writefrequency", po::value<int>(&writefrequency)->default_value(1))
	("savefrequency", po::value<int>(&savefrequency)->default_value(100))
	("sourcefile", po::value<std::string>(&sourcefile)->default_value("conf.00000"))
	("ignore_checksum_errors", po::value<bool>(&ignore_checksum_errors)->default_value(false))
	("print_to_screen", po::value<bool>(&print_to_screen)->default_value(false))
	("host_seed", po::value<uint32_t>(&host_seed)->default_value(4815))
	("initial_prng_state", po::value<std::string>(&initial_prng_state)->default_value(""))

	//gaugefield parameters
	("beta", po::value<double>(&beta)->default_value(4.0))
	("rho", po::value<double>(&rho)->default_value(0.))
	("rho_iter", po::value<int>(&rho_iter)->default_value(0))
	("gaugeact", po::value<std::string>()->default_value("wilson"))

	//heatbath parameters
	("thermalization", po::value<int>(&thermalizationsteps)->default_value(0))
	("heatbathsteps", po::value<int>(&heatbathsteps)->default_value(1000))
	("overrelaxsteps", po::value<int>(&overrelaxsteps)->default_value(1))
	("xi", po::value<int>(&xi)->default_value(1))
	("measure_transportcoefficient_kappa", po::value<bool>(&measure_transportcoefficient_kappa)->default_value(false) )

	//fermionic parameters
	("fermact", po::value<std::string>()->default_value("wilson"))
	("fermact_mp", po::value<std::string>()->default_value("wilson"))
	("kappa", po::value<double>(&kappa)->default_value(0.125))
	("mass", po::value<double>(&mass)->default_value(0.1))
	("mu", po::value<double>(&mu)->default_value(0.006))
	("csw", po::value<double>(&csw)->default_value(0.))
	("kappa_mp", po::value<double>(&kappa_mp)->default_value(0.125))
	("mu_mp", po::value<double>(&mu_mp)->default_value(0.006))
	("csw_mp", po::value<double>(&csw_mp)->default_value(0.))
	("cgmax", po::value<int>(&cgmax)->default_value(1000))
	("cgmax_mp", po::value<int>(&cgmax_mp)->default_value(1000))
	("theta_fermion_spatial", po::value<double>(&theta_fermion_spatial)->default_value(0.))
	("theta_fermion_temporal", po::value<double>(&theta_fermion_temporal)->default_value(0.))
	("chem_pot_re", po::value<double>(&chem_pot_re)->default_value(0.))
	("chem_pot_im", po::value<double>(&chem_pot_im)->default_value(0.))
	("use_eo", po::value<bool>(&use_eo)->default_value(true))
	("solver", po::value<std::string>()->default_value("bicgstab"))
	("solver_mp", po::value<std::string>()->default_value("bicgstab"))
	("use_gauge_only", po::value<bool>(&use_gauge_only)->default_value(false))
	("sourcetype",  po::value<std::string>()->default_value("point"), "Type of source to use for inverter")
	("sourcecontent",  po::value<std::string>()->default_value("one"), "Type of content to use with inverter sources")
	("num_sources", po::value<int>(&num_sources)->default_value(12))
	("source_x", po::value<int>(&source_x)->default_value(0))
	("source_y", po::value<int>(&source_y)->default_value(0))
	("source_z", po::value<int>(&source_z)->default_value(0))
	("source_t", po::value<int>(&source_t)->default_value(0))
	("place_sources_on_host", po::value<bool>(&place_sources_on_host)->default_value(false))
#ifdef _USEDOUBLEPREC_
	("solver_prec", po::value<double>(&solver_prec)->default_value(1e-23))
	("force_prec", po::value<double>(&force_prec)->default_value(1e-12))
#else
	("solver_prec", po::value<double>(&solver_prec)->default_value(1e-16))
	("force_prec", po::value<double>(&force_prec)->default_value(1e-8))
#endif
	("iter_refresh", po::value<int>(&iter_refresh)->default_value(100))
	("iter_refresh_mp", po::value<int>(&iter_refresh_mp)->default_value(100))

	//HMC specific parameters
	("tau", po::value<double>(&tau)->default_value(0.5))
	("reversibility_check", po::value<bool>(&reversibility_check)->default_value(false))
	("integrationsteps0", po::value<int>(&integrationsteps0)->default_value(10))
	("integrationsteps1", po::value<int>(&integrationsteps1)->default_value(10))
	("integrationsteps2", po::value<int>(&integrationsteps2)->default_value(10))
	("hmcsteps", po::value<int>(&hmcsteps)->default_value(10))
	("num_timescales", po::value<int>(&num_timescales)->default_value(1))
	("integrator0", po::value<std::string>()->default_value("leapfrog"))
	("integrator1", po::value<std::string>()->default_value("leapfrog"))
	("integrator2", po::value<std::string>()->default_value("leapfrog"))
	// this is the optimal value...
	("lambda0", po::value<double>(&lambda0)->default_value(0.1931833275037836))
	("lambda1", po::value<double>(&lambda1)->default_value(0.1931833275037836))
	("lambda2", po::value<double>(&lambda2)->default_value(0.1931833275037836))
	
	//RHMC specific parameters
	("md_approx_ord", po::value<int>(&md_approx_ord)->default_value(8))
	("metro_approx_ord", po::value<int>(&metro_approx_ord)->default_value(15))
	("findminmax_max", po::value<int>(&findminmax_max)->default_value(5000))
	("findminmax_iteration_block_size", po::value<int>(&findminmax_iteration_block_size)->default_value(25), "find_minmax will check the residual only every N iterations")
	("findminmax_prec", po::value<double>(&findminmax_prec)->default_value(1.e-3))
	("conservative", po::value<bool>(&conservative)->default_value(false))
	("num_tastes", po::value<int>(&num_tastes)->default_value(2))
	("approx_lower", po::value<double>(&approx_lower)->default_value(1.e-5))
	("approx_upper", po::value<double>(&approx_upper)->default_value(1.))
	("rhmcsteps", po::value<int>(&rhmcsteps)->default_value(10))

	("corr_dir", po::value<int>(&corr_dir)->default_value(3), "Direction for the correlator")

	("use_same_rnd_numbers", po::value<bool>(&use_same_rnd_numbers)->default_value(false), "Use random numbers compatible with a scalar version. SLOW!")
	("profile_solver", po::value<bool>(&profile_solver)->default_value(false))
	("test_ref_val", po::value<double>(&test_ref_value)->default_value(0.))

	("disable-ocl-compiler-opt", po::value<bool>(&ocl_compiler_opt_disabled)->default_value(false), "Disable OpenCL compiler from performing optimizations (adds -cl-disable-opt)")

	("use_merge_kernels_spinor", po::value<bool>(&use_merge_kernels_spinor)->default_value(false), "Use kernel merging for spinor kernels")
	("use_merge_kernels_fermion", po::value<bool>(&use_merge_kernels_fermion)->default_value(false), "Use kernel merging for fermion kernels")
	("use_rec12", po::value<bool>(&use_rec12)->default_value(false), "Use reconstruct 12 compression for SU3 matrices")

	("log-level", po::value<std::string>(&log_level)->default_value("ALL"), "Minimum output log level: ALL TRACE DEBUG INFO WARN ERROR FATAL OFF")

	("read_multiple_configs", po::value<bool>(&read_multiple_configs)->default_value(false), "Read in more than one gaugefield configuration")
	("config_read_start", po::value<int>(&config_read_start)->default_value(0), "Number to begin with when reading in more than one gaugefield configuration")
	("config_read_end", po::value<int>(&config_read_end)->default_value(1), "Number to end with when reading in more than one gaugefield configuration")
	("config_read_incr", po::value<int>(&config_read_incr)->default_value(1), "Increment for gaugefield configuration number when reading in more than one gaugefield configuration")
	("config_number_digits", po::value<int>(&config_number_digits)->default_value(5), "Number of digits to name gaugefield configurations")
	("config_prefix", po::value<std::string>(&config_prefix)->default_value("conf."), "Prefix for gaugefield configuration")
	("config_postfix", po::value<std::string>(&config_postfix)->default_value(""), "Postfix for gaugefield configuration")
	("prng_prefix", po::value<std::string>(&prng_prefix)->default_value("prng."), "Prefix for PRNG configuration")
	("prng_postfix", po::value<std::string>(&prng_postfix)->default_value(""), "Postfix for PRNG configuration")

	//parameters to write out observables
	("gauge_obs_to_single_file", po::value<bool>(&gauge_obs_to_single_file)->default_value(true), "Save gauge observables to one single file")
	("gauge_obs_prefix", po::value<std::string>(&gauge_obs_prefix)->default_value("gaugeObs"), "Prefix for gauge observables file")
	("gauge_obs_postfix", po::value<std::string>(&gauge_obs_postfix)->default_value(".dat"), "Postfix for gauge observables file")
	("ferm_obs_to_single_file", po::value<bool>(&ferm_obs_to_single_file)->default_value(false), "Save fermionic observables to one single file")
	("ferm_obs_corr_prefix", po::value<std::string>(&ferm_obs_corr_prefix)->default_value(""), "Prefix for fermionic observables (correlators) file")
	("ferm_obs_corr_postfix", po::value<std::string>(&ferm_obs_corr_postfix)->default_value("_correlators.dat"), "Postfix for fermionic observables (correlators) file")
	("ferm_obs_pbp_prefix", po::value<std::string>(&ferm_obs_pbp_prefix)->default_value(""), "Prefix for fermionic observables (chiral condensate) file")
	("ferm_obs_pbp_postfix", po::value<std::string>(&ferm_obs_pbp_postfix)->default_value("_pbp.dat"), "Postfix for fermionic observables (chiral condensate) file")
	("hmc_obs_to_single_file", po::value<bool>(&hmc_obs_to_single_file)->default_value(true), "Save hmc observables to one single file")
	("hmc_obs_prefix", po::value<std::string>(&hmc_obs_prefix)->default_value("hmc_output"), "Prefix for hmc observables file")
	("hmc_obs_postfix", po::value<std::string>(&hmc_obs_postfix)->default_value(""), "Postfix for hmc observables file")
	("rhmc_obs_to_single_file", po::value<bool>(&rhmc_obs_to_single_file)->default_value(true), "Save rhmc observables to one single file")
	("rhmc_obs_prefix", po::value<std::string>(&rhmc_obs_prefix)->default_value("rhmc_output"), "Prefix for rhmc observables file")
	("rhmc_obs_postfix", po::value<std::string>(&rhmc_obs_postfix)->default_value(""), "Postfix for rhmc observables file")

	("measure_correlators", po::value<bool>(&measure_correlators)->default_value(true), "Measure fermionic correlators")
	("measure_pbp", po::value<bool>(&measure_pbp)->default_value(false), "Measure chiral condensate")

	("pbp_version",  po::value<std::string>()->default_value("std"), "Version of chiral condensate")

	("cg_iteration_block_size", po::value<int>(&cg_iteration_block_size)->default_value(10), "CG will check the residual only every N iterations")
	("cg_use_async_copy", po::value<bool>(&cg_use_async_copy)->default_value(false), "CG will use residual of iteration N - block_size for termination condition.")
	("cg_minimum_iteration_count", po::value<int>(&cg_minimum_iteration_count)->default_value(0), "CG will perform at least this many itertions. USE ONLY FOR BENCHMARKS!")

	("split_cpu", po::value<bool>(&split_cpu)->default_value(false), "Split the CPU into multiple devices to avoid numa issues. (Requires OpenCL 1.2 at least)");


	po::options_description desc;
	desc.add(cmd_opts).add(config);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).positional(pos_opts).run(), vm);
	if(vm.count("help")) { // see http://stackoverflow.com/questions/5395503/required-and-optional-arguments-using-boost-library-program-options as to why this is done before po::notifiy(vm)
		std::cout << desc << '\n';
		throw Inputparameters::parse_aborted();
	}

	if(vm.count("input-file")) {
		std::string config_file = vm["input-file"].as<std::string>();
		ConfigFileNormalizer normalizer;
		add_option_aliases(&normalizer);
		// add stuff from input file
		std::string normalized_file;
		try {
			normalized_file = normalizer(config_file);
		} catch(std::invalid_argument) {
			std::cout << "Invalid config file " << config_file << std::endl;
			throw Inputparameters::parse_aborted();
		}
		std::istringstream normalized_file_stream(normalized_file);
		po::store(po::parse_config_file(normalized_file_stream, config, false), vm);
	}

	po::notify(vm); // checks whether all required arguments are set

	// handle the enumeration types
	_startcondition = ::get_startcondition(vm["startcondition"].as<std::string>());
	gaugeact = ::get_action(vm["gaugeact"].as<std::string>());
	fermact = ::get_action(vm["fermact"].as<std::string>());
	fermact_mp = ::get_action(vm["fermact_mp"].as<std::string>());
	integrator0 = ::get_integrator(vm["integrator0"].as<std::string>());
	integrator1 = ::get_integrator(vm["integrator1"].as<std::string>());
	integrator2 = ::get_integrator(vm["integrator2"].as<std::string>());
	_solver = ::get_solver(vm["solver"].as<std::string>());
	_solver_mp = ::get_solver(vm["solver_mp"].as<std::string>());
	sourcetype = ::get_sourcetype(vm["sourcetype"].as<std::string>() );
	sourcecontent = ::get_sourcecontent(vm["sourcecontent"].as<std::string>() );
	pbp_version_ = ::get_pbp_version(vm["pbp_version"].as<std::string>() );
}

static Inputparameters::action get_action(std::string s)
{
	boost::algorithm::to_lower(s);
	std::map<std::string, Inputparameters::action> m;
	m["wilson"] = Inputparameters::wilson;
	m["clover"] = Inputparameters::clover;
	m["twistedmass"] = Inputparameters::twistedmass;
	m["tlsym"] = Inputparameters::tlsym;
	m["iwasaki"] = Inputparameters::iwasaki;
	m["dbw2"] = Inputparameters::dbw2;
	m["rooted_stagg"] = Inputparameters::rooted_stagg;

	Inputparameters::action a = m[s];
	if(a) { // map returns 0 if element is not found
		return a;
	} else {
		std::cout << s << " is not a valid action." << std::endl;
		throw Inputparameters::parse_aborted();
	}
}
static Inputparameters::integrator get_integrator(std::string s)
{
	boost::algorithm::to_lower(s);
	std::map<std::string, Inputparameters::integrator> m;
	m["leapfrog"] = Inputparameters::leapfrog;
	m["twomn"] = Inputparameters::twomn;
	m["2mn"] = Inputparameters::twomn;

	Inputparameters::integrator a = m[s];
	if(a) {
		return a;
	} else {
		std::cout << s << " is not a valid integrator." << std::endl;
		throw Inputparameters::parse_aborted();
	}
}
static Inputparameters::startcondition get_startcondition(std::string s)
{
	boost::algorithm::to_lower(s);
	std::map<std::string, Inputparameters::startcondition> m;
	m["cold_start"] = Inputparameters::cold_start;
	m["cold"] = Inputparameters::cold_start;
	m["hot_start"] = Inputparameters::hot_start;
	m["hot"] = Inputparameters::hot_start;
	m["start_from_source"] = Inputparameters::start_from_source;
	m["source"] = Inputparameters::start_from_source;
	m["continue"] = Inputparameters::start_from_source;

	Inputparameters::startcondition a = m[s];
	if(a) {
		return a;
	} else {
		std::cout << s << " is not a valid startcondition." << std::endl;
		throw Inputparameters::parse_aborted();
	}
}
static Inputparameters::solver get_solver(std::string s)
{
	boost::algorithm::to_lower(s);
	std::map<std::string, Inputparameters::solver> m;
	m["cg"] = Inputparameters::cg;
	m["bicgstab"] = Inputparameters::bicgstab;
	m["bicgstab_save"] = Inputparameters::bicgstab_save;
	Inputparameters::solver a = m[s];
	if(a) {
		return a;
	} else {
		std::cout << s << " is not a valid solver." << std::endl;
		throw Inputparameters::parse_aborted();
	}
}
static Inputparameters::sourcetypes get_sourcetype(std::string s)
{
	boost::algorithm::to_lower(s);
	std::map<std::string, Inputparameters::sourcetypes> m;
	m["point"] = Inputparameters::point;
	m["volume"] = Inputparameters::volume;
	m["timeslice"] = Inputparameters::timeslice;
	m["zslice"] = Inputparameters::zslice;

	Inputparameters::sourcetypes a = m[s];
	if(a) { // map returns 0 if element is not found
		return a;
	} else {
		std::cout << s << " is not a valid sourcetype." << std::endl;
		throw Inputparameters::parse_aborted();
	}
}
static Inputparameters::sourcecontents get_sourcecontent(std::string s)
{
	boost::algorithm::to_lower(s);
	std::map<std::string, Inputparameters::sourcecontents> m;
	m["one"] = Inputparameters::one;
	m["z4"] = Inputparameters::z4;
	m["gaussian"] = Inputparameters::gaussian;

	Inputparameters::sourcecontents a = m[s];
	if(a) { // map returns 0 if element is not found
		return a;
	} else {
		std::cout << s << " is not a valid sourcecontent." << std::endl;
		throw Inputparameters::parse_aborted();
	}
}
static Inputparameters::pbp_version get_pbp_version(std::string s)
{
	boost::algorithm::to_lower(s);
	std::map<std::string, Inputparameters::pbp_version> m;
	m["std"] = Inputparameters::std;
	m["tm_one_end_trick"] = Inputparameters::tm_one_end_trick;

	Inputparameters::pbp_version a = m[s];
	if(a) { // map returns 0 if element is not found
		return a;
	} else {
		std::cout << s << " is not a valid pvp_version." << std::endl;
		throw Inputparameters::parse_aborted();
	}
}

static void add_option_aliases(meta::ConfigFileNormalizer * const normalizer)
{
	normalizer->add_alias("NS", "nspace");
	normalizer->add_alias("NT", "ntime");
	normalizer->add_alias("use_evenodd", "use_eo");
	normalizer->add_alias("thermsteps", "thermalization");
	normalizer->add_alias("thermalizationsteps", "thermalization");
	normalizer->add_alias("puregauge", "use_gauge_only");
	normalizer->add_alias("PGT", "use_gauge_only");
	normalizer->add_alias("test_ref_value", "test_ref_val");
	normalizer->add_alias("ThetaT", "theta_fermion_temporal");
	normalizer->add_alias("ThetaS", "theta_fermion_spatial");
}

bool meta::Inputparameters::get_split_cpu() const noexcept
{
	return split_cpu;
}
