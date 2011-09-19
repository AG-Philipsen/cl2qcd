#include "hmc.h"

int main(int argc, char* argv[])
{
  try{

    if(argc != 2) throw Print_Error_Message("Need file name for input parameters",__FILE__,__LINE__);

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
	cl_device_type* devicetypes = new cl_device_type[parameters.get_num_dev()];
	logger.trace() << "init gaugefield" ;
	gaugefield.init(parameters.get_num_dev(), devicetypes, &parameters);
	delete [] devicetypes;
	logger.trace()<< "initial gaugeobservables:";
	gaugefield.print_gaugeobservables(&poly_timer, &plaq_timer);
	size_t rndsize = gaugefield.get_numrndstates();
	init_random_seeds(gaugefield.get_rndarray(), "rand_seeds", rndsize);
	logger.trace() << "Got seeds";
	gaugefield.copy_gaugefield_to_devices();
	gaugefield.copy_rndarray_to_devices();
	logger.trace() << "Moved stuff to device";
	init_timer.add();
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// HMC
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	perform_timer.reset();
	int hmc_iter = parameters.get_hmcsteps();
	int iter;	
	hmc_float acc_rate = 0.;
	int writefreq = parameters.get_writefrequency();
	int savefreq = parameters.get_savefrequency();
	//This is the random-number generator for the metropolis-step
	Random hmc_rnd_gen (parameters.get_host_seed());
	
	logger.trace() << "perform HMC on device... ";
	//main hmc-loop
	for(iter = 0; iter < hmc_iter; iter ++) {
		//generate new random-number for Metropolis step
		hmc_float rnd_number = hmc_rnd_gen.doub();
		gaugefield.perform_hmc_step(0, &parameters, &obs, iter, rnd_number);
		acc_rate += obs.accept;
		if( ( (iter + 1) % writefreq ) == 0 ) {
			gaugefield.print_hmcobservables(obs, iter, gaugeout_name.str());
		}
		if( parameters.get_saveconfigs() == true && ( (iter + 1) % savefreq ) == 0 ) {
			gaugefield.sync_gaugefield();
			gaugefield.save(iter);
		}
	}
	logger.trace() << "HMC done";
	logger.trace() << "Acceptance rate: " << setprecision(1) << percent(acc_rate, hmc_iter) << "%";
	perform_timer.add();
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Final Output
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	total_timer.add();
	general_time_output(&total_timer, &init_timer, &perform_timer, &plaq_timer, &poly_timer);
	//print times from the devices...
// 	logger.trace() << "## Device: HMC";
// 	(gaugefield.get_task_solver())->print_copy_times(totaltime);
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// free variables
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	gaugefield.finalize();

  } //try
  //exceptions from Opencl classes
  catch (Opencl_Error& e) {
    logger.fatal()<<e.what();
    exit(1);
  }
  catch (File_Exception& fe) {
    logger.fatal()<<"Could not open file: "<<fe.get_filename();
    logger.fatal()<<"Aborting.";
    exit(1);
  }
  catch (Print_Error_Message& em) {
    logger.fatal()<<em.what();
    exit(1);
  }

  return 0;
}
