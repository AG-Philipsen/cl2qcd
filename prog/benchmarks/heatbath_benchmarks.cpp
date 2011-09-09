#include "../heatbath.h"

int main(int argc, char* argv[])
{
	
#ifndef _PROFILING_
	logger.fatal() << "_PROFILING_ not defined, cannot perform benchmarks. Aborting...";
	throw Print_Error_Message("_PROFILING_ not defined, cannot perform benchmarks. Aborting...");
#endif

//CP: This should be the same as the normal heatbath-executable
/////////////////////////////////////////////////////////////////////////////////////////

	if(argc != 2) {
		logger.fatal() << "need file name for input parameters";
		throw File_Exception("No file given");
	}

	char* inputfile = argv[1];
	inputparameters parameters;
	parameters.readfile(inputfile);
	parameters.print_info_heatbath(argv[0]);

	//name of file to store gauge observables
	stringstream gaugeout_name;
	gaugeout_name << "gaugeobservables_beta" << parameters.get_beta();

	fstream logfile;
	logfile.open("heatbath.log", std::ios::out | std::ios::app);
	if(logfile.is_open()) {
	  parameters.print_info_heatbath(argv[0], &logfile);
	  logfile.close();
	} else {
	  logger.warn()<<"Could not open heatbath.log";
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Initialization
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	init_timer.reset();
	sourcefileparameters parameters_source;

	Gaugefield_heatbath gaugefield;

	cl_device_type primary_device_type;
	//check whether GPU should be used
	if(parameters.get_use_gpu() == true) {
		primary_device_type = CL_DEVICE_TYPE_GPU;
	} else {
		primary_device_type = CL_DEVICE_TYPE_CPU;
	}
	gaugefield.init(1, primary_device_type, &parameters);
	logger.trace() << "Got gaugefield";
	gaugefield.print_gaugeobservables(0);

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
		gaugefield.perform_tasks(parameters.get_overrelaxsteps());
		gaugefield.synchronize(0);
		gaugefield.print_gaugeobservables_from_task(i, 0, gaugeout_name.str());
	}
	logger.trace() << "heatbath-benchmarking done";
	perform_timer.add();

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Final Output
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//TODO: remove gaugeobservables-file, this is not really needed
	total_timer.add();
	general_time_output(&total_timer, &init_timer, &perform_timer, (gaugefield.get_task_heatbath())->get_copy_to(), gaugefield.get_task_heatbath()->get_copy_on(), &plaq_timer, &poly_timer);

	//CP: this is just a fist version and will go into an own file later
	stringstream profiling_out;
	profiling_out << argv[0] << "_profiling_data";

	fstream prof_file;
	prof_file.open(profiling_out.str(), std::ios::out | std::ios::app);
	if(prof_file.is_open()) {
	  parameters.print_info_heatbath(argv[0], &prof_file);
	  prof_file.close();
	} else {
	  logger.warn()<<"Could not open " << profiling_out;
	}
	gaugefield.print_profiling(profiling_out.str());
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// free variables
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	gaugefield.finalize();

	return 0;
}
