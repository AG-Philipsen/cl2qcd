#include "hmc.h"

int main(int argc, char* argv[]) {

  char* progname = argv[0];
  print_hello(progname);

  char* inputfile = argv[1];
  inputparameters parameters;
  parameters.readfile(inputfile);
  print_info(&parameters);

  benchmark_id = argv[2];
  
  stringstream gaugeout_name;
  gaugeout_name<<"gaugeobservables_beta"<<parameters.get_beta();
  
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Initialization
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

  sourcefileparameters parameters_source;
  hmc_gaugefield * gaugefield;
  gaugefield = (hmc_gaugefield*) malloc(sizeof(hmc_gaugefield));
  hmc_rndarray rndarray;

  init_gaugefield(gaugefield,&parameters,&inittime);
  init_random_seeds(rnd, rndarray, &inittime);

	simple_correlator(gaugefield, parameters.get_kappa(), parameters.get_mu(), parameters.get_theta_fermion(), parameters.get_chem_pot_re(), parameters.get_chem_pot_im(), 1000);
	
#ifdef _USEGPU_
  opencl device(CL_DEVICE_TYPE_GPU, local_work_size, global_work_size, &inittime);
#else
  opencl device(CL_DEVICE_TYPE_CPU, local_work_size, global_work_size, &inittime);
#endif

  cout << "initial values of observables:\n\t" ;
  print_gaugeobservables(gaugefield, &polytime, &plaqtime);

  device.copy_gaugefield_to_device(gaugefield, &copytime);
  device.copy_rndarray_to_device(rndarray, &copytime);

#ifdef _TESTING_
  device.testing(gaugefield);
#endif
	usetimer noop;
	device.init_solver_variables(&parameters, local_work_size, global_work_size, &noop);
	device.testing_spinor(&parameters, local_work_size, global_work_size);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Heatbath
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////

	int nsteps = parameters.get_heatbathsteps();
	int overrelaxsteps = parameters.get_overrelaxsteps();
	cout<<"perform "<<nsteps<<" heatbath steps on OpenCL device..."<<endl;
	for(int i = 0; i<nsteps; i++){
		device.run_heatbath(parameters.get_beta(), local_work_size, global_work_size, &updatetime);
		for (int j = 0; j < overrelaxsteps; j++)
			device.run_overrelax(parameters.get_beta(), local_work_size, global_work_size, &overrelaxtime);
		if( ( (i+1)%parameters.get_writefrequency() ) == 0 ) {
			device.gaugeobservables(local_work_size, global_work_size, &plaq, &tplaq, &splaq, &pol, &plaqtime, &polytime);
			print_gaugeobservables(plaq, tplaq, splaq, pol, i, gaugeout_name.str());
		}
		if( parameters.get_saveconfigs()==TRUE && ( (i+1)%parameters.get_savefrequency() ) == 0 ) {
			device.get_gaugefield_from_device(gaugefield, &copytime);
			save_gaugefield(gaugefield, &parameters, i);
			print_gaugeobservables(gaugefield, &plaqtime, &polytime, i, gaugeout_name.str());
		}
	}

	device.get_gaugefield_from_device(gaugefield, &copytime);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Final Output
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  totaltime.add();
  save_gaugefield(gaugefield, &parameters, nsteps);  
  time_output(&totaltime, &inittime, &polytime, &plaqtime, &updatetime, &overrelaxtime, &copytime);

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // free variables
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  free(gaugefield);
  device.finalize();
  
  return HMC_SUCCESS;
}
