#include "../inverter.h"

int main(int argc, char* argv[])
{
	try {

		if(argc != 2) {
			logger.fatal() << "need file name for input parameters";
			throw File_Exception("No file given");
		}

		char* inputfile = argv[1];
		inputparameters parameters;
		parameters.readfile(inputfile);
		parameters.print_info_inverter(argv[0]);

		ofstream ofile;
		ofile.open("inverter.log");
		if(ofile.is_open()) {
			parameters.print_info_inverter(argv[0], &ofile);
			ofile.close();
		} else {
			logger.warn() << "Could not log file for inverter.";
		}

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Initialization
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		init_timer.reset();
		Gaugefield_inverter gaugefield;

		//use 2 devices: one for solver, one for correlator
		int numtasks = 2;
		if(parameters.get_num_dev() != 2 )
			logger.warn() << "Only 1 device demanded by input file. All calculations performed on primary device.";

		cl_device_type primary_device;
		switch ( parameters.get_use_gpu() ) {
			case true :
				primary_device = CL_DEVICE_TYPE_GPU;
				break;
			case false :
				primary_device = CL_DEVICE_TYPE_CPU;
				break;
		}

		logger.trace() << "Init gaugefield" ;
		gaugefield.init(numtasks, primary_device, &parameters);


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
			/** @todo Insert functions here */
		}
		logger.trace() << "inverter-benchmarking done" ;
		perform_timer.add();

		//////////////////////////////////////////////////////////////////////////////////////////////////////////////
		// Final Output
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////

		total_timer.add();
		uint64_t totaltime = total_timer.getTime();
		general_time_output(&total_timer, &init_timer, &perform_timer, &plaq_timer, &poly_timer);
		//print times from the devices...
		logger.info() << "## Device: Solver";
		(gaugefield.get_task_solver())->print_copy_times(totaltime);
		logger.info() << "## Device: Correlator";
		(gaugefield.get_task_correlator())->print_copy_times(totaltime);
		
		//CP: this is just a fist version and will go into an own file later
		stringstream profiling_out;
		profiling_out << argv[0] << "_profiling_data";

		fstream prof_file;
		prof_file.open(profiling_out.str(), std::ios::out | std::ios::app);
		if(prof_file.is_open()) {
			parameters.print_info_inverter(argv[0], &prof_file);
			prof_file.close();
		} else {
			logger.warn() << "Could not open " << profiling_out;
		}
		gaugefield.print_profiling(profiling_out.str());

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
