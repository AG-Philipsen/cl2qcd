#include "inverter.h"

#include "physics/lattices/gaugefield.hpp"
#include "physics/lattices/spinorfield.hpp"
#include "physics/sources.hpp"
#include "physics/algorithms/flavour_doublet.hpp"
#include "physics/algorithms/inversion.hpp"

#include "meta/util.hpp"

void perform_measurements(hardware::System& system, physics::lattices::Gaugefield& gf, physics::PRNG& prng, const std::string config_name);

int main(int argc, const char* argv[])
{
	using namespace physics::lattices;
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
			logger.warn() << "Could not open log file for inverter.";
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Initialization
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		init_timer.reset();
		hardware::System system(parameters);
		physics::PRNG prng(system);

		init_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// inverter
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		perform_timer.reset();
		logger.info() << "Perform inversion(s) on device.." ;

		if(parameters.get_read_multiple_configs()) {
			const int iter_start = (parameters.get_read_multiple_configs()) ? parameters.get_config_read_start() : 0;
			// the "+1" is necessary so that all config are read
			const int iter_end   = (parameters.get_read_multiple_configs()) ? parameters.get_config_read_end() + 1   : 1;
			const int iter_incr =  (parameters.get_read_multiple_configs()) ? parameters.get_config_read_incr()  : 1;
			int iter = 0;

			//main loop
			for(iter = iter_start; iter < iter_end; iter += iter_incr) {
				std::string config_name = meta::create_configuration_name(parameters, iter);
				logger.info() << config_name;
				logger.info() << "Measure fermionic observables on configuration: " << config_name;
				Gaugefield gaugefield(system, prng, config_name) ;
				perform_measurements(system, gaugefield, prng, config_name);
			}
		} else {
			std::string config_name = parameters.get_sourcefile();
			logger.info() << config_name;
			logger.info() << "Measure fermionic observables on configuration: " << config_name;
			Gaugefield gaugefield(system, prng);
			perform_measurements(system, gaugefield, prng, config_name);
		}
		logger.trace() << "Inversion done" ;
		perform_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Final Output
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		total_timer.add();
		general_time_output(&total_timer, &init_timer, &perform_timer, &plaq_timer, &poly_timer);

		if(parameters.get_profile_solver() ) {
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
		}

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
	} catch (Invalid_Parameters& es) {
		logger.fatal() << es.what();
		exit(1);
	}

	return 0;

}

void perform_measurements(hardware::System& system, physics::lattices::Gaugefield& gaugefield, physics::PRNG& prng, const std::string config_name)
{
	using namespace physics;
	using namespace physics::lattices;
	using namespace physics::algorithms;

	auto parameters = system.get_inputparameters();
	if(parameters.get_print_to_screen() ) {
		print_gaugeobservables(gaugefield, gaugefield.get_parameters_source().trajectorynr_source);
	}
	if(parameters.get_measure_correlators() ) {
		// for the correlator calculation, all sources are needed on the device
		const std::vector<Spinorfield*> sources = create_swappable_sources(system, prng, parameters.get_num_sources());
		const std::vector<Spinorfield*> result = create_swappable_spinorfields(system, sources.size(), parameters.get_place_sources_on_host());

		swap_out(sources);
		swap_out(result);

		perform_inversion(&result, &gaugefield, sources, system);

		logger.info() << "Finished inversion. Starting measurements.";

		swap_in(sources);
		swap_in(result);

		//get name for file to which correlators are to be stored
		std::string corr_fn = meta::get_ferm_obs_corr_file_name(parameters, config_name);
		flavour_doublet_correlators(result, sources, corr_fn, parameters);
		release_spinorfields(result);
		release_spinorfields(sources);
	}
	if(parameters.get_measure_pbp() ) {
		//get name for file to which pbp is to be stored
		std::string pbp_fn = meta::get_ferm_obs_pbp_file_name(parameters, config_name);
		// the chiral condensate needs only one source at a time
		for(int i_sources = 0; i_sources < parameters.get_num_sources(); i_sources ++) {
			auto sources = create_sources(system, prng, 1);
			auto result = create_spinorfields(system, sources.size());

			perform_inversion(&result, &gaugefield, sources, system);

			flavour_doublet_chiral_condensate(result, sources, pbp_fn, gaugefield.get_parameters_source().trajectorynr_source, system);
			release_spinorfields(result);
			release_spinorfields(sources);
		}
	}
}
