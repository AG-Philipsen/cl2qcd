/** @file
 * Input file handling implementation
 *
 * Copyright (c) 2012 Matthias Bach <bach@compeng.uni-frankfurt.de>
 * Copyright (c) 2012, 2014 Christopher Pinke <pinke@compeng.uni-frankfurt.de>
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

bool Inputparameters::get_use_merge_kernels_fermion() const noexcept
{
  return use_merge_kernels_fermion;
}
bool Inputparameters::get_use_merge_kernels_spinor() const noexcept
{
  return use_merge_kernels_spinor;
}

Inputparameters::sourcetypes Inputparameters::get_sourcetype() const noexcept
{
  return sourcetype;
}
Inputparameters::sourcecontents Inputparameters::get_sourcecontent() const noexcept
{
  return sourcecontent;
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
	
	po::options_description options_gaugefield("Gaugefield options");
	options_gaugefield.add_options()
		("ignore_checksum_errors", po::value<bool>(&ignore_checksum_errors)->default_value(false))
		("beta", po::value<double>(&beta)->default_value(4.0))
		("use_smearing", po::value<bool>(&use_smearing)->default_value(false))
		("rho", po::value<double>(&rho)->default_value(0.))
		("rho_iter", po::value<int>(&rho_iter)->default_value(0))
		("gaugeact", po::value<std::string>()->default_value("wilson"));
	
	po::options_description options_heatbath("Heatbath options");
	options_heatbath.add_options()
		//todo: this is also needed in the HMC!
		("thermalization", po::value<int>(&thermalizationsteps)->default_value(0))
		("heatbathsteps", po::value<int>(&heatbathsteps)->default_value(1000))
		("overrelaxsteps", po::value<int>(&overrelaxsteps)->default_value(1))
		//todo: is this used?
		("xi", po::value<int>(&xi)->default_value(1));

	po::options_description options_fermion("Fermion options");
	options_fermion.add_options()
		("fermact", po::value<std::string>()->default_value("wilson"))
		("fermact_mp", po::value<std::string>()->default_value("wilson"))
		//todo: change this default value!
		("kappa", po::value<double>(&kappa)->default_value(0.125))
		("mass", po::value<double>(&mass)->default_value(0.1))
		("mu", po::value<double>(&mu)->default_value(0.006))
		("csw", po::value<double>(&csw)->default_value(0.))
		("kappa_mp", po::value<double>(&kappa_mp)->default_value(0.125))
		("mu_mp", po::value<double>(&mu_mp)->default_value(0.006))
		("csw_mp", po::value<double>(&csw_mp)->default_value(0.))
		("theta_fermion_spatial", po::value<double>(&theta_fermion_spatial)->default_value(0.))
		("theta_fermion_temporal", po::value<double>(&theta_fermion_temporal)->default_value(0.))
		("chem_pot_re", po::value<double>(&chem_pot_re)->default_value(0.))
		("chem_pot_im", po::value<double>(&chem_pot_im)->default_value(0.))
		("use_chem_pot_re", po::value<bool>(&use_chem_pot_re)->default_value(false))
		("use_chem_pot_im", po::value<bool>(&use_chem_pot_im)->default_value(false))
		("use_eo", po::value<bool>(&use_eo)->default_value(true))
		("use_merge_kernels_spinor", po::value<bool>(&use_merge_kernels_spinor)->default_value(false), "Use kernel merging for spinor kernels")
		("use_merge_kernels_fermion", po::value<bool>(&use_merge_kernels_fermion)->default_value(false), "Use kernel merging for fermion kernels");


	po::options_description options_source("Source options");
	options_source.add_options()
		("sourcetype",  po::value<std::string>()->default_value("point"), "Type of source to use for inverter")
		("sourcecontent",  po::value<std::string>()->default_value("one"), "Type of content to use with inverter sources")
		("num_sources", po::value<int>(&num_sources)->default_value(12))
		("source_x", po::value<int>(&source_x)->default_value(0))
		("source_y", po::value<int>(&source_y)->default_value(0))
		("source_z", po::value<int>(&source_z)->default_value(0))
		("source_t", po::value<int>(&source_t)->default_value(0))
		("place_sources_on_host", po::value<bool>(&place_sources_on_host)->default_value(false));

	po::options_description options_solver("Solver options");
	options_solver.add_options()
		("solver", po::value<std::string>()->default_value("bicgstab"))
		("solver_mp", po::value<std::string>()->default_value("bicgstab"))
		("cgmax", po::value<int>(&cgmax)->default_value(1000))
		("cgmax_mp", po::value<int>(&cgmax_mp)->default_value(1000))
#ifdef _USEDOUBLEPREC_
		("solver_prec", po::value<double>(&solver_prec)->default_value(1e-23))
		("force_prec", po::value<double>(&force_prec)->default_value(1e-12))
#else
		("solver_prec", po::value<double>(&solver_prec)->default_value(1e-16))
		("force_prec", po::value<double>(&force_prec)->default_value(1e-8))
#endif
		("iter_refresh", po::value<int>(&iter_refresh)->default_value(100))
		("iter_refresh_mp", po::value<int>(&iter_refresh_mp)->default_value(100))
		//this is not used. Remove!
		("profile_solver", po::value<bool>(&profile_solver)->default_value(false))
		("cg_iteration_block_size", po::value<int>(&cg_iteration_block_size)->default_value(10), "CG will check the residual only every N iterations")
		("cg_use_async_copy", po::value<bool>(&cg_use_async_copy)->default_value(false), "CG will use residual of iteration N - block_size for termination condition.")
		("cg_minimum_iteration_count", po::value<int>(&cg_minimum_iteration_count)->default_value(0), "CG will perform at least this many itertions. USE ONLY FOR BENCHMARKS!");


	
	po::options_description desc;
	desc.add(cmd_opts)
		.add(ParametersConfig::getOptions())
		.add(ParametersIo::getOptions())
		.add(options_gaugefield).add(options_heatbath).add(options_fermion).add(options_source).add(options_solver)
		.add(ParametersObs::getOptions())
		.add(ParametersHmc::getOptions())
		.add(ParametersRhmc::getOptions())
		.add(ParametersTest::getOptions());
	
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
		//todo: do not pass the "help" option...
		po::store(po::parse_config_file(normalized_file_stream, desc, false), vm);
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
	m["z2"] = Inputparameters::z2;

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
	normalizer->add_alias("test_ref_value2", "test_ref_val2");
	normalizer->add_alias("ThetaT", "theta_fermion_temporal");
	normalizer->add_alias("ThetaS", "theta_fermion_spatial");
}
