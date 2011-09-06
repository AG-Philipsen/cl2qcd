#include "tk_kappa_hybrid.h"

int main(int argc, char* argv[])
{
  try {

    if(argc != 2) throw Print_Error_Message("Need file name for input parameters",__FILE__,__LINE__);

	char* inputfile = argv[1];
	inputparameters parameters;
	parameters.readfile(inputfile);
	parameters.print_info_tkkappa(argv[0]);

	//name of file to store gauge observables
	stringstream gaugeout_name;
	gaugeout_name << "gaugeobservables_beta" << parameters.get_beta();


	fstream logfile;	
	logfile.open("tk_kappa_hybrid.log", std::ios::out | std::ios::app);
	if(logfile.is_open()) {
	  parameters.print_info_tkkappa(argv[0],&logfile);
	  logfile.close();
	} else {
	  throw File_Exception("tk_kappa_hybrid.log");
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Initialization
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////


	Gaugefield_heatbath_kappa gaugefield;
	int numtasks = 2;

	// this is the device type for the heatbath
	cl_device_type primary_device_type = CL_DEVICE_TYPE_GPU;

        gaugefield.init(numtasks, primary_device_type, &parameters);

	gaugefield.print_gaugeobservables(0);
	gaugefield.print_gaugeobservables(0, gaugeout_name.str());
	for(int iter = 0; iter < parameters.get_heatbathsteps() / parameters.get_writefrequency(); iter++) {
	  gaugefield.perform_tasks(parameters.get_writefrequency(), parameters.get_overrelaxsteps());
	  gaugefield.synchronize(0);
	  gaugefield.print_gaugeobservables(iter);
	  //	  gaugefield.print_gaugeobservables_from_task(iter,0);
	  //	  gaugefield.print_gaugeobservables_from_task(iter,1);
	  gaugefield.print_gaugeobservables(iter,gaugeout_name.str());
	  //	  gaugefield.print_kappa(iter,"kappa_clover.dat");
	}

	gaugefield.save("conf.save");

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
