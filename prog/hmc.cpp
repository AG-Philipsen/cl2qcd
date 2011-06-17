#include "hmc.h"

int main(int argc, char* argv[])
{

	if(argc != 2) {
		logger.fatal() << "need file name for input parameters";
		return HMC_FILEERROR;
	}

	char* progname = argv[0];
	print_hello(progname);

	char* inputfile = argv[1];
	inputparameters parameters;
	parameters.readfile(inputfile);
	print_info(&parameters,&cout);

	//init file to store gauge observables, print initial information
	stringstream gaugeout_name;
	gaugeout_name << "gaugeobservables_beta" << parameters.get_beta();
	fstream gaugeout;

	gaugeout.open(gaugeout_name.str().c_str(), std::ios::out | std::ios::app);
	if(!gaugeout.is_open()) exit(HMC_FILEERROR);
	print_info(&parameters, &gaugeout);
	gaugeout.close();	
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Initialization
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	sourcefileparameters parameters_source;

	Gaugefield_hmc gaugefield;
	hmc_rndarray rndarray;
	cl_device_type devicetypes[1];

#ifdef _USEGPU_
	devicetypes[0] = CL_DEVICE_TYPE_GPU;
#else
	devicetypes[0] = CL_DEVICE_TYPE_CPU;
#endif
	cerr << "init gaugefield" << endl;
	gaugefield.init(1, devicetypes, &parameters, &inittime);
	logger.trace() << "Got gaugefield";
	
	/** @todo this needs to be implemented using structs.. */
	//cerr << "print initial gaugeobservables..." << endl;
	//	gaugefield.print_gaugeobservables(&polytime, &plaqtime);

	int err = init_random_seeds(rndarray, "rand_seeds", &inittime);
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
	//main hmc-loop
	logger.trace() << "start main HMC loop with " << hmc_iter << " iterations: " ;
	for(iter = 0; iter < hmc_iter; iter ++) {
		//generate new random-number for Metropolis step
		rnd_number = hmc_rnd_gen.doub();
		gaugefield.perform_hmc_step(&parameters, iter, rnd_number,  &copytimer,&singletimer,&Mtimer,&scalarprodtimer,&latimer,&dslashtimer,&Mdiagtimer,&solvertimer);
		if( ( (iter + 1) % writefreq ) == 0 ) {
 			gaugefield.print_gaugeobservables_from_devices(&plaqtime, &polytime, iter, gaugeout_name.str());
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
