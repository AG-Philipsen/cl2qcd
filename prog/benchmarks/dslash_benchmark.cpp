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

		//one needs 2 task here
		int numtasks = 2;
		if(parameters.get_device_count() != numtasks )
			logger.warn() << "Exactly 2 devices demanded by benchmark executable. All calculations performed on primary device.";

		cl_device_type primary_device = parameters.get_use_gpu() ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU;

		logger.trace() << "Init gaugefield" ;
		gaugefield.init(numtasks, primary_device);

		logger.info() << "Gaugeobservables:";
		gaugefield.print_gaugeobservables(0);

		//init needed buffers again
		//these are: 2 eoprec spinorfield, 1 gaugefield
		const hardware::buffers::Spinor sf1(meta::get_eoprec_spinorfieldsize(parameters), gaugefield.get_task_solver()->get_device());
		const hardware::buffers::Spinor sf2(meta::get_eoprec_spinorfieldsize(parameters), gaugefield.get_task_solver()->get_device());

		init_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// dslash-benchmark
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		perform_timer.reset();

		int hmc_iter = parameters.get_hmcsteps();
		int iter;
		auto gf_buffer = gaugefield.get_task_solver()->get_device()->get_gaugefield_code()->get_gaugefield();

		logger.trace() << "Perform " << hmc_iter << "of dslash benchmarking (EVEN + ODD) for each step";
		for(iter = 0; iter < hmc_iter; iter ++) {
			gaugefield.get_task_solver()->dslash_eo_device(&sf1, &sf2, gf_buffer, EVEN);
			gaugefield.get_task_solver()->dslash_eo_device(&sf1, &sf2, gf_buffer, ODD);
		}
		logger.trace() << "dslash benchmarking done" ;
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

		gaugefield.get_task_solver()->print_profiling(profiling_out, 0) ;

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
