#include "tk_kappa.h"

int main(int argc, char* argv[])
{

	if(argc != 2) {
		cerr << "need file name for input parameters" << endl;
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
	print_info(&parameters,&gaugeout);
	gaugeout.close();
	
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
	if(ntherm > 0) gaugefield.heatbath(local_work_size, global_work_size, ntherm, &updatetime);

	int nsteps = parameters.get_heatbathsteps();
	int overrelaxsteps = parameters.get_overrelaxsteps();
	int writefreq = parameters.get_writefrequency();
	int savefreq = parameters.get_savefrequency();
	
	hmc_float kappa_karsch_val;
	hmc_float kappa_clover_val;
	ofstream kappa_karsch_out;
	kappa_karsch_out.open ("kappa_karsch.dat");
	ofstream kappa_clover_out;
	kappa_clover_out.open ("kappa_clover.dat");

	for(int i = 0; i < nsteps; i++) {
		gaugefield.heatbath(local_work_size, global_work_size, &updatetime);
		for(int j = 0; j < overrelaxsteps; j++) gaugefield.overrelax(local_work_size, global_work_size, &overrelaxtime);
		if( ( (i + 1) % writefreq ) == 0 ) {
			gaugefield.print_gaugeobservables_from_devices(local_work_size, global_work_size, &plaqtime, &polytime, i, gaugeout_name.str());
		}
		if( parameters.get_saveconfigs() == TRUE && ( (i + 1) % savefreq ) == 0 ) {
			gaugefield.sync_gaugefield(&copytime);
			gaugefield.save(i);
		}
	//Add a measurement frequency
	gaugefield.sync_gaugefield(&copytime);
	hmc_error err;
	err = gaugefield.kappa_karsch ();
	err = gaugefield.kappa_clover ();

	kappa_karsch_out << gaugefield.get_kappa_karsch() <<endl;
	kappa_clover_out << gaugefield.get_kappa_clover() <<endl;
	}
	
	kappa_karsch_out.close();
	kappa_clover_out.close();
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
