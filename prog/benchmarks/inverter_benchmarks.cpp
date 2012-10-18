#include "../inverter.h"

#include "../meta/util.hpp"

#include <boost/program_options.hpp>
namespace po = boost::program_options;

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
			logger.warn() << "Could not log file for inverter.";
		}

		//get name for file to which correlators are to be stored
		stringstream corr_fn;
		switch ( parameters.get_startcondition() ) {
			case meta::Inputparameters::start_from_source :
				corr_fn << parameters.get_sourcefile() << "_correlators.dat" ;
				break;
			case meta::Inputparameters::hot_start :
				corr_fn << "conf.hot_correlators.dat" ;
				break;
			case meta::Inputparameters::cold_start :
				corr_fn << "conf.cold_correlators.dat" ;
				break;
		}

		if(parameters.get_profile_solver() == false) {
			logger.warn() << "solver times will not be measured!";
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Initialization
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		init_timer.reset();
		hardware::System system(parameters, true);
		Gaugefield_inverter gaugefield(&system);

		//one needs 2 tasks here since the correlator-module produces the sources...
		int numtasks = 2;
		if(parameters.get_device_count() != 2 )
			logger.warn() << "Only 1 device demanded by benchmark executable. All calculations performed on primary device.";

		cl_device_type primary_device = parameters.get_use_gpu() ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU;

		logger.trace() << "Init gaugefield" ;
		gaugefield.init(numtasks, primary_device);


		logger.info() << "Gaugeobservables:";
		gaugefield.print_gaugeobservables(0);
		init_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// inverter-benchmarks
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		perform_timer.reset();

		int hmc_iter = parameters.get_hmcsteps();
		int iter;

		logger.trace() << "Perform " << hmc_iter << "of benchmarking";
		for(iter = 0; iter < hmc_iter; iter ++) {
			//CP: these are esssentially the same actions as the "normal" inverter performs...
			logger.info() << "Perform inversion on device.." ;

			gaugefield.create_sources();
			gaugefield.perform_inversion(&solver_timer);

			//flavour_doublet_correlators does a sync at the beginning
			gaugefield.flavour_doublet_correlators(corr_fn.str());

			logger.trace() << "Inversion done" ;

		}
		logger.trace() << "inverter-benchmarking done" ;

		if(parameters.get_profile_solver() == true) {
			logger.info() << "Inverter took " << solver_timer.getTime() << " ms";
		}

		perform_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Final Output
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		total_timer.add();
		uint64_t totaltime = total_timer.getTime();
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
		gaugefield.print_profiling(profiling_out);

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
	}
}
