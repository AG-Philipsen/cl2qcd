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
static meta::action get_action(std::string);
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
	
	po::options_description options_heatbath("Heatbath options");
	options_heatbath.add_options()
		//todo: this is also needed in the HMC!
		("thermalization", po::value<int>(&thermalizationsteps)->default_value(0))
		("heatbathsteps", po::value<int>(&heatbathsteps)->default_value(1000))
		("overrelaxsteps", po::value<int>(&overrelaxsteps)->default_value(1))
		//todo: is this used?
		("xi", po::value<int>(&xi)->default_value(1));

	po::options_description desc;
	desc.add(cmd_opts)
		.add(ParametersConfig::getOptions())
		.add(ParametersIo::getOptions())
		.add(ParametersGauge::getOptions())
		.add(options_heatbath)
		.add(ParametersFermion::getOptions())
		.add(ParametersSolver::getOptions())
		.add(ParametersSources::getOptions())
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

static meta::action get_action(std::string s)
{
	boost::algorithm::to_lower(s);
	std::map<std::string, meta::action> m;
	m["wilson"] = meta::action::wilson;
	m["clover"] = meta::action::clover;
	m["twistedmass"] = meta::action::twistedmass;
	m["tlsym"] = meta::action::tlsym;
	m["iwasaki"] = meta::action::iwasaki;
	m["dbw2"] = meta::action::dbw2;
	m["rooted_stagg"] = meta::action::rooted_stagg;

	meta::action a = m[s];
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
