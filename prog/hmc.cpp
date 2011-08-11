#include "hmc.h"

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

	sourcefileparameters parameters_source;
	hmc_observables obs;
	
	Gaugefield_hmc gaugefield;
	hmc_rndarray rndarray;
	cl_device_type devicetypes[parameters.get_num_dev()];
	gaugefield.init_devicetypes_array(devicetypes, &parameters);

	logger.trace() << "init gaugefield" ;
	gaugefield.init(parameters.get_num_dev(), devicetypes, &parameters);
	logger.trace()<< "initial gaugeobservables:";
	gaugefield.print_gaugeobservables(&polytime, &plaqtime);

	int err = init_random_seeds(rndarray, "rand_seeds");
	if(err) return err;
	logger.trace() << "Got seeds";

	gaugefield.copy_gaugefield_to_devices(&copytimer);
	gaugefield.copy_rndarray_to_devices(rndarray, &copytime);
	
	logger.trace() << "Moved stuff to device";
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// HMC
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	int hmc_iter = parameters.get_hmcsteps();
	int iter;	
	int writefreq = parameters.get_writefrequency();
	int savefreq = parameters.get_savefrequency();
	hmc_float rnd_number;
	/** @todo CP: give seed meaningful value, perhaps read it in from parameters */
	int seed = 89343894;
	Random hmc_rnd_gen (seed);
	
	logger.trace() << "perform HMC on device... ";
	logger.trace() << "start main HMC loop with " << hmc_iter << " iterations: " ;
	//main hmc-loop
	for(iter = 0; iter < hmc_iter; iter ++) {
		//generate new random-number for Metropolis step
		rnd_number = hmc_rnd_gen.doub();
		gaugefield.perform_hmc_step(0, &parameters, &obs, iter, rnd_number, gaugeout_name.str(),  &copytimer,&singletimer,&Mtimer,&scalarprodtimer,&latimer,&dslashtimer,&Mdiagtimer,&solvertimer);
		if( ( (iter + 1) % writefreq ) == 0 ) {
			gaugefield.print_hmcobservables(obs, iter, gaugeout_name.str());
		}
		if( parameters.get_saveconfigs() == TRUE && ( (iter + 1) % savefreq ) == 0 ) {
			gaugefield.sync_gaugefield(&copytime);
			gaugefield.save(iter);
		}
	}
	logger.trace() << "HMC done";
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// free variables
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	err = gaugefield.finalize();
	if (err!= HMC_SUCCESS) 
		logger.fatal() << "error in finalizing HMC";
	return HMC_SUCCESS;
}
