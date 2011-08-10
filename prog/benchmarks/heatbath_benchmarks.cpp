//CP: this is done in a first attempt
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

	sourcefileparameters parameters_source;

	Gaugefield_heatbath gaugefield;
	hmc_rndarray rndarray;
	cl_device_type devicetypes[parameters.get_num_dev()];

	if(parameters.get_num_dev() == 1){
	#ifdef _USEGPU_
		devicetypes[0] = CL_DEVICE_TYPE_GPU;
	#else
		devicetypes[0] = CL_DEVICE_TYPE_CPU;
	#endif
	}
	else if(parameters.get_num_dev() == 2){
		devicetypes[0] = CL_DEVICE_TYPE_GPU;
		devicetypes[1] = CL_DEVICE_TYPE_CPU;
	}
	else{
		logger.fatal() << "Number of devices too big, aborting..." ;
		return HMC_STDERR;
	}

	inittime.reset();
	gaugefield.init(1, devicetypes, &parameters);
	inittime.add();
	
	logger.trace() << "Got gaugefield";

	int err = init_random_seeds(rndarray, "rand_seeds", &inittime);
	if(err) return err;

	logger.trace() << "Got seeds";

	//first output, if you like it...
	logger.trace() << "OpenCL initialisaton time:\t" << inittime.getTime() << " [mus]";
	gaugefield.print_gaugeobservables(&polytime,&plaqtime);

	gaugefield.copy_gaugefield_to_devices(&copytime);
	gaugefield.copy_rndarray_to_devices(rndarray, &copytime);

	logger.trace() << "Moved stuff to device";

/////////////////////////////////////////////////////////////////////////////////////////	
//CP: Now it differs from the normal heatbath-executable

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Heatbath-benchmark
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int nsteps = parameters.get_heatbathsteps();
	
	logger.info() << "Perform " << nsteps << "of benchmarking";

	for(int i = 0; i < nsteps; i++) {
		gaugefield.heatbath();
		gaugefield.overrelax();
		gaugefield.print_gaugeobservables_from_devices(&plaqtime, &polytime, i, gaugeout_name.str(), 0);
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Final Output
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	//TODO: remove gaugeobservables-file
	
	//CP: this is just a fist version and will go into an own file later
	stringstream profiling_out;
	profiling_out << "profiling";

	fstream prof_file;
	prof_file.open(profiling_out.str(), std::ios::out | std::ios::app);
	if(prof_file.is_open()) {
	  parameters.print_info_heatbath(argv[0], &prof_file);
	  prof_file.close();
	} else {
	  logger.warn()<<"Could not open " << profiling_out;
	}
	gaugefield.get_devices()[0].print_profiling(profiling_out.str());

	totaltime.add();
	time_output_heatbath(&totaltime, &inittime, &polytime, &plaqtime, &updatetime, &overrelaxtime, &copytime);


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// free variables
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	gaugefield.finalize();

	return HMC_SUCCESS;
}
