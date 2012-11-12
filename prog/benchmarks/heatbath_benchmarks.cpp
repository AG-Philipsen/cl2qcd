#include "../general_header.h"

#include "../meta/util.hpp"
#include "../physics/prng.hpp"
#include "../physics/lattices/gaugefield.hpp"
#include "../physics/algorithms/heatbath.hpp"

#include "../meta/util.hpp"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc, const char* argv[])
{
//CP: This should be the same as the normal heatbath-executable
/////////////////////////////////////////////////////////////////////////////////////////
	using namespace physics;
	using namespace physics::lattices;
	using physics::algorithms::heatbath;

	meta::Inputparameters parameters(argc, argv);
	switchLogLevel(parameters.get_log_level());

	meta::print_info_heatbath(argv[0], parameters);

	//name of file to store gauge observables
	stringstream gaugeout_name;
	gaugeout_name << "gaugeobservables_beta" << parameters.get_beta();

	fstream logfile;
	logfile.open("heatbath.log", std::ios::out | std::ios::app);
	if(logfile.is_open()) {
		meta::print_info_heatbath(argv[0], &logfile, parameters);
		logfile.close();
	} else {
		logger.warn() << "Could not open heatbath.log";
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Initialization
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	init_timer.reset();

	hardware::System system(parameters, true);

	PRNG prng(system);
	Gaugefield gaugefield(system, prng);

	print_gaugeobservables(gaugefield, 0);

	init_timer.add();

/////////////////////////////////////////////////////////////////////////////////////////
//CP: Now it differs from the normal heatbath-executable

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Heatbath-benchmark
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	perform_timer.reset();
	int nsteps = parameters.get_heatbathsteps();

	logger.info() << "Perform " << nsteps << "of benchmarking";

	for(int i = 0; i < nsteps; i++) {
		heatbath(gaugefield, prng, parameters.get_overrelaxsteps());
		print_gaugeobservables(gaugefield, i, gaugeout_name.str());
	}
	logger.trace() << "heatbath-benchmarking done";
	perform_timer.add();

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Final Output
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//TODO: remove gaugeobservables-file, this is not really needed
	total_timer.add();
	general_time_output(&total_timer, &init_timer, &perform_timer, &plaq_timer, &poly_timer);

	//CP: this is just a fist version and will go into an own file later
	string profiling_out;
	profiling_out = string(argv[0]) + string("_profiling_data");

	fstream prof_file;
	prof_file.open(profiling_out.c_str(), std::ios::out | std::ios::app);
	if(prof_file.is_open()) {
		meta::print_info_heatbath(argv[0], &prof_file, parameters);
		prof_file.close();
	} else {
		logger.warn() << "Could not open " << profiling_out;
	}

	print_profiling(system, profiling_out);

	return 0;
}
