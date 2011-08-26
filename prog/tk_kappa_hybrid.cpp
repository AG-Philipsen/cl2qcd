#include "tk_kappa_hybrid.h"

int main(int argc, char* argv[])
{

	if(argc != 2) {
	  logger.warn() << "need file name for input parameters" ;
	  return HMC_FILEERROR;
	}

	char* inputfile = argv[1];
	inputparameters parameters;
	parameters.readfile(inputfile);
	parameters.print_info_tkkappa(argv[0]);

	//name of file to store observables
	stringstream gaugeout_name;
	gaugeout_name << "gaugeobservables_beta" << parameters.get_beta();
       	stringstream kappa_karsch_out;
	kappa_karsch_out<<"kappa_karsch.dat";
       	stringstream kappa_clover_out;
	kappa_karsch_out<<"kappa_clover.dat";

	fstream logfile;	
	logfile.open("tk_kappa_hybrid.log", std::ios::out | std::ios::app);
	if(logfile.is_open()) {
	  parameters.print_info_tkkappa(argv[0],&logfile);
	  logfile.close();
	} else {
	  logger.warn() << "Could not open tk_kappa.log";
	  exit(HMC_FILEERROR);
	}
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Initialization
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	sourcefileparameters parameters_source;

	Gaugefield_k_hybrid gaugefield;
	cl_device_type devicetypes[2];

	devicetypes[0] = CL_DEVICE_TYPE_GPU;
	devicetypes[1] = CL_DEVICE_TYPE_CPU;

	int num_ocl_devices[2] = {1,1};

	gaugefield.init(num_ocl_devices, 2, devicetypes, &parameters);

	gaugefield.copy_gaugefield_to_devices(&copytime);


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Heatbath
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////


	int ntherm = parameters.get_thermalizationsteps();
	if(ntherm > 0) gaugefield.heatbath(ntherm, &updatetime);

	int nsteps = parameters.get_heatbathsteps();
	int overrelaxsteps = parameters.get_overrelaxsteps();
	int writefreq = parameters.get_writefrequency();	//LZ Note: Use writefreq to determine frequency for kappa measurements 
	int savefreq = parameters.get_savefrequency();
	

	logger.trace() << "Start heatbath and measurement of TK kappa" ;
	
	usetimer timer_karsch;
	usetimer timer_clover;
	hmc_float time_karsch = 0.;
	hmc_float time_clover = 0.;
	
	
	for(int i = 0; i < nsteps/writefreq; i++) {

	  {
	    
	    for(int n=0; n<writefreq; n++) {
	      gaugefield.heatbath(&updatetime);
	      for(int j = 0; j < overrelaxsteps; j++)
		gaugefield.overrelax(&overrelaxtime);
	    }
	    cout<<"print"<<endl;
      	    gaugefield.print_gaugeobservables_from_devices(&plaqtime, &polytime, i, gaugeout_name.str());

	    hmc_error err;
	    err = gaugefield.kappa_karsch_gpu (&timer_karsch);
	    err |= gaugefield.kappa_clover_gpu (&timer_clover);
	    gaugefield.print_TK(kappa_karsch_out.str(),kappa_clover_out.str());
	    
	    time_karsch += timer_karsch.getTime();
	    time_clover += timer_clover.getTime();
	  
	  } //end OMP sections
	  gaugefield.sync_gaugefield(&copytime);	    
	  if( parameters.get_saveconfigs() && ( (i + 1) % savefreq ) == 0 ) 
	    gaugefield.save(i);

 	}
	

	logger.info() <<"Measurement TK kappa_karsch: " << time_karsch/1000000. / hmc_float (nsteps) << " s" ;
	logger.info() <<"Measurement TK kappa_clover: " << time_clover/1000000. / hmc_float (nsteps) << " s" ;
	
	gaugefield.sync_gaugefield(&copytime);
	gaugefield.save(nsteps);
	
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Final Output
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	totaltime.add();
	// FIXME time_output_heatbath(&totaltime, &inittime, &polytime, &plaqtime, &updatetime, &overrelaxtime, &copytime);


	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// free variables
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	gaugefield.finalize();

	return HMC_SUCCESS;
}
