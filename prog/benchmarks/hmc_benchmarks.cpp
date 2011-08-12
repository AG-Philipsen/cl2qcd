#include "../hmc.h"

int main(int argc, char* argv[])
{

	if(argc != 2) {
		logger.fatal() << "need file name for input parameters";
		return HMC_FILEERROR;
	}

	char* inputfile = argv[1];
	inputparameters parameters;
	parameters.readfile(inputfile);
	parameters.print_info_hmc(argv[0]);

	//name of file to store gauge observables, print initial information
	/** @todo think about what is a senseful filename*/
	stringstream gaugeout_name;
	gaugeout_name << "hmc_output";

	fstream logfile;
	logfile.open("hmc.log", std::ios::out | std::ios::app);
	if(logfile.is_open()) {
	  parameters.print_info_hmc(argv[0],&logfile);
	  logfile.close();	
	} else {
	  logger.warn() << "Could not open hmc.log";
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Initialization
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	init_timer.reset();
	sourcefileparameters parameters_source;
	hmc_observables obs;
	
	Gaugefield_hmc gaugefield;
	hmc_rndarray rndarray;
	cl_device_type devicetypes[parameters.get_num_dev()];
	gaugefield.init_devicetypes_array(devicetypes, &parameters);
	logger.trace() << "init gaugefield" ;
	gaugefield.init(parameters.get_num_dev(), devicetypes, &parameters);
	logger.trace()<< "initial gaugeobservables:";
	gaugefield.print_gaugeobservables(&poly_timer, &plaq_timer);
	int err = init_random_seeds(rndarray, "rand_seeds");
	if(err) return err;
	logger.trace() << "Got seeds";
	gaugefield.copy_gaugefield_to_devices(&copy_to_from_dev_timer);
	gaugefield.copy_rndarray_to_devices(rndarray, &copy_to_from_dev_timer);
	logger.trace() << "Moved stuff to device";
	init_timer.add();
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// HMC-Benchmarks
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	perform_timer.reset();
	int hmc_iter = parameters.get_hmcsteps();
	int iter;	
	//This is the random-number generator for the metropolis-step
// 	Random hmc_rnd_gen (parameters.get_host_seed());
	
	logger.trace() << "Perform " << hmc_iter << "of benchmarking";
	for(iter = 0; iter < hmc_iter; iter ++) {
		/** @todo Insert functions here */
	}
	logger.trace() << "HMC-benchmarking done";
	perform_timer.add();
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Final Output
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
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

	err = gaugefield.finalize();
	if (err!= HMC_SUCCESS) 
		logger.fatal() << "error in finalizing " << argv[0];
	return HMC_SUCCESS;
}
