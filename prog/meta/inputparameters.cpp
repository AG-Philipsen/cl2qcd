/** @file
 * Input file handling implementation
 */

#include "inputparameters.h"

#include <map>
#include <iostream>
#include <fstream>
#include <boost/program_options.hpp>
#include <boost/algorithm/string.hpp>
namespace po = boost::program_options;

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
bool Inputparameters::get_saveconfigs() const noexcept
{
	return saveconfigs;
}
int Inputparameters::get_writefrequency() const noexcept
{
	return writefrequency;
}
int Inputparameters::get_savefrequency() const noexcept
{
	return savefrequency;
}
std::string Inputparameters::get_sourcefilenumber() const noexcept
{
	return sourcefilenumber;
}
bool Inputparameters::get_print_to_screen() const noexcept
{
	return print_to_screen;
}
//This is obvious!!!
uint64_t Inputparameters::get_host_seed() const noexcept
{
	return host_seed;
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
double Inputparameters::get_mu() const noexcept
{
	return mu;
}
double Inputparameters::get_csw() const noexcept
{
	return csw;
}
int Inputparameters::get_iter0() const noexcept
{
	return iter0;
}
int Inputparameters::get_iter1() const noexcept
{
	return iter1;
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
int Inputparameters::get_iter0_mp() const noexcept
{
	return iter0_mp;
}
int Inputparameters::get_iter1_mp() const noexcept
{
	return iter1_mp;
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
//at the moment, only 2 solvers are implemented..
bool Inputparameters::get_use_cg() const noexcept
{
	return use_cg;
}
bool Inputparameters::get_use_cg_mp() const noexcept
{
	return use_cg_mp;
}
bool Inputparameters::get_use_bicgstab_save() const noexcept
{
	return use_bicgstab_save;
}
bool Inputparameters::get_use_bicgstab_save_mp() const noexcept
{
	return use_bicgstab_save_mp;
}
bool Inputparameters::get_use_pointsource() const noexcept
{
	return use_pointsource;
}
bool Inputparameters::get_use_gauge_only() const noexcept
{
	return use_gauge_only;
}
int Inputparameters::get_num_sources() const noexcept
{
	return num_sources;
}
int Inputparameters::get_pointsource_x() const noexcept
{
	return pointsource_x;
}
int Inputparameters::get_pointsource_y() const noexcept
{
	return pointsource_y;
}
int Inputparameters::get_pointsource_z() const noexcept
{
	return pointsource_z;
}
int Inputparameters::get_pointsource_t() const noexcept
{
	return pointsource_t;
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

//HMC specific parameters
double Inputparameters::get_tau() const noexcept
{
	return tau;
}
bool Inputparameters::get_reversibility_check() const noexcept
{
	return reversibility_check;
}
int Inputparameters::get_integrationsteps0() const noexcept
{
	return integrationsteps0;
}
int Inputparameters::get_integrationsteps1() const noexcept
{
	return integrationsteps1;
}
int Inputparameters::get_hmcsteps() const noexcept
{
	return hmcsteps;
}
int Inputparameters::get_num_timescales() const noexcept
{
	return num_timescales;
}
Inputparameters::integrator Inputparameters::get_integrator0() const noexcept
{
	return integrator0;
}
Inputparameters::integrator Inputparameters::get_integrator1() const noexcept
{
	return integrator1;
}
Inputparameters::integrator Inputparameters::get_integrator2() const noexcept
{
	return integrator2;
}
//this is the optimal value...
double Inputparameters::get_lambda0() const noexcept
{
	return lambda0;
}
double Inputparameters::get_lambda1() const noexcept
{
	return lambda1;
}
double Inputparameters::get_lambda2() const noexcept
{
	return lambda2;
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

Inputparameters::Inputparameters(int argc, const char** argv)
{
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

	("use_aniso", po::value<bool>(&use_aniso)->default_value(false))
	("use_chem_pot_re", po::value<bool>(&use_chem_pot_re)->default_value(false))
	("use_chem_pot_im", po::value<bool>(&use_chem_pot_im)->default_value(false))
	("use_smearing", po::value<bool>(&use_smearing)->default_value(false))
	("use_mp", po::value<bool>(&use_mp)->default_value(false))
	("nspace", po::value<int>(&nspace)->default_value(4))
	("ntime", po::value<int>(&ntime)->default_value(8))

	("startcondition", po::value<std::string>()->default_value("start_from_source"))
	("saveconfigs", po::value<bool>(&saveconfigs)->default_value(false))
	("writefrequency", po::value<int>(&writefrequency)->default_value(1))
	("savefrequency", po::value<int>(&savefrequency)->default_value(100))
	("sourcefile", po::value<std::string>(&sourcefilenumber)->default_value("00000"))
	("print_to_screen", po::value<bool>(&print_to_screen)->default_value(false))
	//This is obvious!!!
	("host_seed", po::value<uint64_t>(&host_seed)->default_value(4815))

	//gaugefield parameters
	("beta", po::value<double>(&beta)->default_value(4.0))
	("rho", po::value<double>(&rho)->default_value(0.))
	("rho_iter", po::value<int>(&rho_iter)->default_value(0))
	("gaugeact", po::value<std::string>()->default_value("wilson"))

	//heatbath parameters
	("thermalizationsteps", po::value<int>(&thermalizationsteps)->default_value(0))
	("heatbathsteps", po::value<int>(&heatbathsteps)->default_value(1000))
	("overrelaxsteps", po::value<int>(&overrelaxsteps)->default_value(1))
	("xi", po::value<int>(&xi)->default_value(1))

	//fermionic parameters
	("fermact", po::value<std::string>()->default_value("wilson"))
	("fermact_mp", po::value<std::string>()->default_value("wilson"))
	("kappa", po::value<double>(&kappa)->default_value(0.125))
	("mu", po::value<double>(&mu)->default_value(0.006))
	("csw", po::value<double>(&csw)->default_value(0.))
	("iter0", po::value<int>(&iter0)->default_value(0))
	("iter1", po::value<int>(&iter1)->default_value(0))
	("kappa_mp", po::value<double>(&kappa_mp)->default_value(0.125))
	("mu_mp", po::value<double>(&mu_mp)->default_value(0.006))
	("csw_mp", po::value<double>(&csw_mp)->default_value(0.))
	("iter0_mp", po::value<int>(&iter0_mp)->default_value(0))
	("iter1_mp", po::value<int>(&iter1_mp)->default_value(0))
	("cgmax", po::value<int>(&cgmax)->default_value(1000))
	("cgmax_mp", po::value<int>(&cgmax_mp)->default_value(1000))
	("theta_fermion_spatial", po::value<double>(&theta_fermion_spatial)->default_value(0.))
	("theta_fermion_temporal", po::value<double>(&theta_fermion_temporal)->default_value(0.))
	("chem_pot_re", po::value<double>(&chem_pot_re)->default_value(0.))
	("chem_pot_im", po::value<double>(&chem_pot_im)->default_value(0.))
	("use_eo", po::value<bool>(&use_eo)->default_value(true))
	// at the moment, only 2 solvers are implemented..
	("use_cg", po::value<bool>(&use_cg)->default_value(false)) // TODO maybe better have solver = ...
	("use_cg_mp", po::value<bool>(&use_cg_mp)->default_value(false))
	("use_bicgstab_save", po::value<bool>(&use_bicgstab_save)->default_value(false))
	("use_bicgstab_save_mp", po::value<bool>(&use_bicgstab_save_mp)->default_value(false))
	("use_pointsource", po::value<bool>(&use_pointsource)->default_value(true))
	("use_gauge_only", po::value<bool>(&use_gauge_only)->default_value(false))
	("num_sources", po::value<int>(&num_sources)->default_value(12))
	("pointsource_x", po::value<int>(&pointsource_x)->default_value(0))
	("pointsource_y", po::value<int>(&pointsource_y)->default_value(0))
	("pointsource_z", po::value<int>(&pointsource_z)->default_value(0))
	("pointsource_t", po::value<int>(&pointsource_t)->default_value(0))
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
	("hmcsteps", po::value<int>(&hmcsteps)->default_value(10))
	("num_timescales", po::value<int>(&num_timescales)->default_value(1))
	("integrator0", po::value<std::string>()->default_value("leapfrog"))
	("integrator1", po::value<std::string>()->default_value("leapfrog"))
	("integrator2", po::value<std::string>()->default_value("leapfrog"))
	// this is the optimal value...
	("lambda0", po::value<double>(&lambda0)->default_value(0.1931833275037836))
	("lambda1", po::value<double>(&lambda1)->default_value(0.1931833275037836))
	("lambda2", po::value<double>(&lambda2)->default_value(0.1931833275037836))

	("corr_dir", po::value<int>(&corr_dir)->default_value(3), "Direction for the correlator")

	("use_same_rnd_numbers", po::value<bool>(&use_same_rnd_numbers)->default_value(false), "Use random numbers compatible with a scalar version. SLOW!")
	("profile_solver", po::value<bool>(&profile_solver)->default_value(false));

	po::options_description desc;
	desc.add(cmd_opts).add(config);

	po::variables_map vm;
	po::store(po::command_line_parser(argc, argv).options(desc).positional(pos_opts).run(), vm);
	if(vm.count("help")) { // see http://stackoverflow.com/questions/5395503/required-and-optional-arguments-using-boost-library-program-options as to why this is done before po::notifiy(vm)
		std::cout << desc << '\n';
		throw Inputparameters::parse_aborted();
	}

	if(vm.count("input-file")) {
		// add stuff from input file
		const std::string config_file_name = vm["input-file"].as<std::string>();
		std::ifstream config_file(config_file_name.c_str());
		if(!config_file.is_open()) {
			std::cout << "Failed to open file " << config_file_name << std::endl;
			throw Inputparameters::parse_aborted();
		}

		// TODO add custom parser to handle other writings for the given parameters

		po::store(po::parse_config_file(config_file, config, false), vm);
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
	m["hot_start"] = Inputparameters::hot_start;
	m["start_from_source"] = Inputparameters::start_from_source;
	m["continue"] = Inputparameters::start_from_source;

	Inputparameters::startcondition a = m[s];
	if(a) {
		return a;
	} else {
		std::cout << s << " is not a valid startcondition." << std::endl;
		throw Inputparameters::parse_aborted();
	}
}
