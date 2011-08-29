#include "../heatbath.h"

int main(int argc, char* argv[])
{
	
#ifndef _PROFILING_
	logger.fatal() << "_PROFILING_ not defined, cannot perform benchmarks. Aborting...";
	exit (HMC_STDERR);
#endif

//CP: This should be the same as the normal heatbath-executable
/////////////////////////////////////////////////////////////////////////////////////////

	if(argc != 2) {
		logger.fatal() << "need file name for input parameters";
		return HMC_FILEERROR;
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

	cl_int err;

	init_timer.reset();
	sourcefileparameters parameters_source;

	Gaugefield_heatbath gaugefield;
	cl_device_type * devicetypes = new cl_device_type[parameters.get_num_dev()];
	gaugefield.init_devicetypes_array(devicetypes, &parameters);

	gaugefield.init(1, devicetypes, &parameters);
	logger.trace() << "Got gaugefield";
	gaugefield.print_gaugeobservables(&poly_timer,&plaq_timer);
	gaugefield.copy_gaugefield_to_devices(&copy_to_from_dev_timer);
	logger.trace() << "Moved stuff to device";
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
		gaugefield.heatbath();
		gaugefield.overrelax();
		gaugefield.print_gaugeobservables_from_devices(i, gaugeout_name.str(), parameters.get_print_to_screen());
	}
	logger.trace() << "heatbath-benchmarking done";
	perform_timer.add();

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Final Output
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	//TODO: remove gaugeobservables-file, this is not really needed
	total_timer.add();
	general_time_output(&total_timer, &init_timer, &perform_timer, &copy_to_from_dev_timer, &copy_on_dev_timer, &plaq_timer, &poly_timer);

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
	gaugefield.get_devices()[0].print_profiling(profiling_out.str());
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// free variables
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	delete[] devicetypes;

	err = gaugefield.finalize();
	if (err!= HMC_SUCCESS) 
		logger.fatal() << "error in finalizing " << argv[0];
	return HMC_SUCCESS;
	
}
