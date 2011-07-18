#include "heatbath.h"

int main(int argc, char* argv[])
{

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


	if(parameters.get_perform_heatbath() != 1){
		logger.fatal() << "perform_heatbath is set to a value other than 1. Aborting...";
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Initialization
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	sourcefileparameters parameters_source;

	Gaugefield gaugefield;
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

	gaugefield.init(1, devicetypes, &parameters, &inittime);

	logger.trace() << "Got gaugefield";

	int err = init_random_seeds(rndarray, "rand_seeds", &inittime);
	if(err) return err;

	logger.trace() << "Got seeds";

	//first output, if you like it...
	//  cout << endl << "OpenCL initialisaton time:\t" << inittime.getTime() << " [mus]" << endl;
	//  gaugefield.print_gaugeobservables(&polytime,&plaqtime);

	gaugefield.copy_gaugefield_to_devices(&copytime);
	gaugefield.copy_rndarray_to_devices(rndarray, &copytime);

	logger.trace() << "Moved stuff to device";

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Heatbath
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	cout<< "Start heatbath" <<endl;

	int ntherm = parameters.get_thermalizationsteps();
	if(ntherm > 0) gaugefield.heatbath(ntherm, &updatetime);

	int nsteps = parameters.get_heatbathsteps();
	int overrelaxsteps = parameters.get_overrelaxsteps();
	int writefreq = parameters.get_writefrequency();
	int savefreq = parameters.get_savefrequency();

	logger.info() << "Start heatbath";

	for(int i = 0; i < nsteps; i++) {
	  
		gaugefield.heatbath(&updatetime);
		for(int j = 0; j < overrelaxsteps; j++) gaugefield.overrelax(&overrelaxtime);
		if( ( (i + 1) % writefreq ) == 0 ) {
		  gaugefield.print_gaugeobservables_from_devices(&plaqtime, &polytime, i, gaugeout_name.str());
		}
		if( parameters.get_saveconfigs() == TRUE && ( (i + 1) % savefreq ) == 0 ) {
			gaugefield.sync_gaugefield(&copytime);
			gaugefield.save(i);
		}
	}

  	gaugefield.sync_gaugefield(&copytime);
 	gaugefield.save(nsteps);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Final Output
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	totaltime.add();
	time_output_heatbath(&totaltime, &inittime, &polytime, &plaqtime, &updatetime, &overrelaxtime, &copytime);


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// free variables
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	gaugefield.finalize();

	return HMC_SUCCESS;
}
