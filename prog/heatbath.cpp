#include "heatbath.h"

int main(int argc, char* argv[])
{

  try {

    if(argc != 2) throw Print_Error_Message("Need file name for input parameters",__FILE__,__LINE__);

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
	cl_device_type* devicetypes = new cl_device_type[parameters.get_num_dev()];
	gaugefield.init(1, devicetypes, &parameters);
	delete [] devicetypes;

	logger.trace() << "Got gaugefield";
	size_t rndsize = gaugefield.get_numrndstates();
	init_random_seeds(gaugefield.get_rndarray(), "rand_seeds", rndsize);
	logger.trace() << "Got seeds";
	gaugefield.print_gaugeobservables(&poly_timer,&plaq_timer);
	gaugefield.copy_gaugefield_to_devices();
	gaugefield.copy_rndarray_to_devices();
	logger.trace() << "Moved stuff to device";
	init_timer.add();

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Heatbath
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	perform_timer.reset();
	logger.trace() << "Start thermalization" ;

	int ntherm = parameters.get_thermalizationsteps();
	if(ntherm > 0) gaugefield.heatbath(ntherm);

	int nsteps = parameters.get_heatbathsteps();
	int overrelaxsteps = parameters.get_overrelaxsteps();
	int writefreq = parameters.get_writefrequency();
	int savefreq = parameters.get_savefrequency();

	logger.info() << "Start heatbath";

	for(int i = 0; i < nsteps; i++) {
		gaugefield.heatbath();
		for(int j = 0; j < overrelaxsteps; j++) gaugefield.overrelax();
		if( ( (i + 1) % writefreq ) == 0 ) {
		  gaugefield.print_gaugeobservables_from_devices(i, gaugeout_name.str(), parameters.get_print_to_screen());
		}
		if( parameters.get_saveconfigs()==true && ( (i + 1) % savefreq ) == 0 ) {
			gaugefield.sync_gaugefield();
			gaugefield.save(i);
		}
	}

  gaugefield.sync_gaugefield();
 	gaugefield.save(nsteps);
	logger.trace() << "heatbath done";
	perform_timer.add();

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Final Output
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	total_timer.add();
	general_time_output(&total_timer, &init_timer, &perform_timer, gaugefield.get_devices_heatbath()[0].get_copy_to(), gaugefield.get_devices_heatbath()[0].get_copy_on(), &plaq_timer, &poly_timer);

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
