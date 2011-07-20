#include "tk_kappa.h"

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

	//name of file to store gauge observables
	stringstream gaugeout_name;
	gaugeout_name << "gaugeobservables_beta" << parameters.get_beta();


	fstream logfile;	
	logfile.open("tk_kappa.log", std::ios::out | std::ios::app);
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

	Gaugefield_k gaugefield;
	hmc_rndarray rndarray;
	cl_device_type devicetypes[1];

#ifdef _USEGPU_
	devicetypes[0] = CL_DEVICE_TYPE_GPU;
#else
	devicetypes[0] = CL_DEVICE_TYPE_CPU;
#endif

	gaugefield.init(1, devicetypes, &parameters, &inittime);
	int err = init_random_seeds(rndarray, "rand_seeds", &inittime);
	if(err) return err;

	//first output, if you like it...
	//  cout << endl << "OpenCL initialisaton time:\t" << inittime.getTime() << " [mus]" << endl;
	//  gaugefield.print_gaugeobservables(&polytime,&plaqtime);

	gaugefield.copy_gaugefield_to_devices(&copytime);
	gaugefield.copy_rndarray_to_devices(rndarray, &copytime);



	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Heatbath
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int ntherm = parameters.get_thermalizationsteps();
	if(ntherm > 0) gaugefield.heatbath(ntherm, &updatetime);

	int nsteps = parameters.get_heatbathsteps();
	int overrelaxsteps = parameters.get_overrelaxsteps();
	int writefreq = parameters.get_writefrequency();
	int savefreq = parameters.get_savefrequency();
	
	ofstream kappa_karsch_out;
	kappa_karsch_out.open ("kappa_karsch.dat");
	kappa_karsch_out.precision(15);
	ofstream kappa_clover_out;
	kappa_clover_out.open ("kappa_clover.dat");
	kappa_clover_out.precision(15);
// 	ofstream q_plaq_out;
// 	q_plaq_out.open ("Q_plaquette.dat");
// 	q_plaq_out.precision(15);

	logger.trace() << "Start heatbath and measurement of TK kappa" ;
	
	usetimer timer_karsch;
	usetimer timer_clover;
	hmc_float time_karsch = 0.;
	hmc_float time_clover = 0.;
	
	
	for(int i = 0; i < nsteps; i++) {
		gaugefield.heatbath(&updatetime);
		for(int j = 0; j < overrelaxsteps; j++)
		  gaugefield.overrelax(&overrelaxtime);
		if( ( (i + 1) % writefreq ) == 0 ) {
			gaugefield.print_gaugeobservables_from_devices(&plaqtime, &polytime, i, gaugeout_name.str());
		}
		if( parameters.get_saveconfigs() == TRUE && ( (i + 1) % savefreq ) == 0 ) {
			gaugefield.sync_gaugefield(&copytime);
			gaugefield.save(i);
		}
	//Add a measurement frequency
	

	//GPU
	hmc_error err;
	err = gaugefield.kappa_karsch_gpu (&timer_karsch);
	err = gaugefield.kappa_clover_gpu (&timer_clover);
	
	//CPU
	gaugefield.sync_gaugefield(&copytime);
//  	err = gaugefield.kappa_karsch ();
//  	err = gaugefield.kappa_clover ();

// 	hmc_float qplaq = gaugefield.Q_plaquette();
// 	q_plaq_out << qplaq <<endl;

	kappa_karsch_out << gaugefield.get_kappa_karsch()  <<endl;
	kappa_clover_out << gaugefield.get_kappa_clover()  <<endl;
	
	time_karsch += timer_karsch.getTime();
	time_clover += timer_clover.getTime();
	
	}
		
	cout.precision(4);
	cout <<"Measurement TK kappa_karsch: " << time_karsch/1000000. / hmc_float (nsteps) << " s"<< endl;
	cout <<"Measurement TK kappa_clover: " << time_clover/1000000. / hmc_float (nsteps) << " s" <<endl;
	
	kappa_karsch_out.close();
	kappa_clover_out.close();
// 	q_plaq_out.close();
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
