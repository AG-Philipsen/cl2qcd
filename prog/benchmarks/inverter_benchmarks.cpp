#include "../inverter.h"

#include "../meta/util.hpp"
#include "../physics/algorithms/inversion.hpp"
#include "../physics/algorithms/flavour_doublet.hpp"
#include "../physics/sources.hpp"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

int main(int argc, const char* argv[])
{
	using namespace physics::lattices;
	using namespace physics::algorithms;
	using namespace physics;

	try {
		meta::Inputparameters parameters(argc, argv);
		switchLogLevel(parameters.get_log_level());

		meta::print_info_inverter(argv[0], parameters);

		ofstream ofile;
		ofile.open("inverter.log");
		if(ofile.is_open()) {
			meta::print_info_inverter(argv[0], &ofile, parameters);
			ofile.close();
		} else {
			logger.warn() << "Could not log file for inverter.";
		}

		//get name for file to which correlators are to be stored
		string corr_fn;
		switch ( parameters.get_startcondition() ) {
			case meta::Inputparameters::start_from_source :
				corr_fn = parameters.get_sourcefile() + "_correlators.dat" ;
				break;
			case meta::Inputparameters::hot_start :
				corr_fn = "conf.hot_correlators.dat" ;
				break;
			case meta::Inputparameters::cold_start :
				corr_fn = "conf.cold_correlators.dat" ;
				break;
		}

		logger.warn() << "solver times will not be measured!";

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Initialization
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		init_timer.reset();
		hardware::System system(parameters, true);
		physics::PRNG prng(system);
		Gaugefield gaugefield(system, prng);


		logger.info() << "Gaugeobservables:";
		print_gaugeobservables(gaugefield, 0);


		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// inverter-benchmarks
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		perform_timer.reset();

		int hmc_iter = parameters.get_hmcsteps();
		int iter;

		{
			ofstream corr_file(corr_fn.c_str(), ios_base::app);
			if(!corr_file.is_open()) {
				throw File_Exception(corr_fn);
			}

			logger.trace() << "Perform " << hmc_iter << "of benchmarking";
			for(iter = 0; iter < hmc_iter; iter ++) {
				//CP: these are esssentially the same actions as the "normal" inverter performs...
				logger.info() << "Perform inversion on device.." ;

				const std::vector<const Spinorfield*> sources = create_sources(system, prng);
				const std::vector<const Spinorfield*> result = create_spinorfields(system, sources.size());
				flavour_doublet_correlators(result, sources, corr_file, parameters);

				logger.trace() << "Inversion done" ;
				release_spinorfields(result);
				release_spinorfields(sources);
			}
			logger.trace() << "inverter-benchmarking done" ;
		}

		perform_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Final Output
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		total_timer.add();
		general_time_output(&total_timer, &init_timer, &perform_timer, &plaq_timer, &poly_timer);

		//CP: this is just a fist version and will go into an own file later
		string profiling_out;
		profiling_out = string(argv[0]) + string("_profiling_data");

		fstream prof_file;
		prof_file.open(profiling_out.c_str(), std::ios::out | std::ios::app);
		if(prof_file.is_open()) {
			meta::print_info_inverter(argv[0], &prof_file, parameters);
			prof_file.close();
		} else {
			logger.warn() << "Could not open " << profiling_out;
		}
		print_solver_profiling(profiling_out);
		print_profiling(system, profiling_out);

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// free variables
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	} //try
	//exceptions from Opencl classes
	catch (Opencl_Error& e) {
		logger.fatal() << e.what();
		exit(1);
	} catch (File_Exception& fe) {
		logger.fatal() << "Could not open file: " << fe.get_filename();
		logger.fatal() << "Aborting.";
		exit(1);
	} catch (Print_Error_Message& em) {
		logger.fatal() << em.what();
		exit(1);
	}
}
