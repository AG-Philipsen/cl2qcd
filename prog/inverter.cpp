#include "inverter.h"

#include "meta/util.hpp"

int main(int argc, const char* argv[])
{
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
		Gaugefield_inverter gaugefield(&system);

		//use 2 devices: one for solver, one for correlator
		int numtasks = 2;
		if(parameters.get_device_count() != 2 )
			logger.warn() << "Only 1 device demanded by input file. All calculations performed on primary device.";

		cl_device_type primary_device = parameters.get_use_gpu() ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU;

		//check if correlator-device is a GPU and in that case exit because the kernels are not meant to be executed there
		if ( parameters.get_use_gpu() == false && parameters.get_device_count() == 2) {
			throw Print_Error_Message("GPU cannot be used for correlator-calculation.", __FILE__, __LINE__);
		}

		logger.trace() << "Init gaugefield" ;
		gaugefield.init(numtasks, primary_device);

		init_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// inverter
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		perform_timer.reset();
		/** @todo usage of solver_timer has to be checked. No output yet */
		usetimer solver_timer;

		logger.info() << "Perform inversion(s) on device.." ;

		if(parameters.get_read_multiple_configs()) {
			int iter_end = parameters.get_config_read_end();
			int iter_start = parameters.get_config_read_start();
			int iter_incr = parameters.get_config_read_incr();
			int iter = 0;

			//main loop
			for(iter = iter_start; iter < iter_end; iter += iter_incr) {
			  std::string config_name = meta::create_configuration_name(parameters, iter);
				logger.info() << "Measure fermionic observables on configuration: " << config_name;
				gaugefield.init_gaugefield(config_name.c_str());
				gaugefield.synchronize(0);
				if(parameters.get_print_to_screen() ) {
					gaugefield.print_gaugeobservables(iter);
				}
				gaugefield.create_sources();
				gaugefield.perform_inversion(&solver_timer);

				//get name for file to which correlators are to be stored
				std::string corr_fn = meta::get_ferm_obs_file_name(parameters, config_name);

				//flavour_doublet_correlators does a sync at the beginning
				gaugefield.flavour_doublet_correlators(corr_fn);
			}
		} else {
			logger.info() << "Gaugeobservables:";
			gaugefield.print_gaugeobservables(0);

			gaugefield.create_sources();
			gaugefield.perform_inversion(&solver_timer);

			//get name for file to which correlators are to be stored
			std::string corr_fn = meta::get_ferm_obs_file_name(parameters, "");

			//flavour_doublet_correlators does a sync at the beginning
			gaugefield.flavour_doublet_correlators(corr_fn);
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

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// free variables
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		gaugefield.finalize();

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
