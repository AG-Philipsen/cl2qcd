#include "hmc.h"

int main(int argc, char* argv[]) {

  char* progname = argv[0];
  print_hello(progname);

  char* inputfile = argv[1];
  inputparameters parameters;
  parameters.readfile(inputfile);
  print_info(&parameters);

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

#ifndef _PERFORM_BENCHMARKS_

#ifdef _FERMIONS_
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Fermions
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	cout<<"calculate simple_correlator on host..."<<endl;
	simple_correlator(gaugefield, parameters.get_kappa(), parameters.get_mu(), parameters.get_theta_fermion(), parameters.get_chem_pot_re(), parameters.get_chem_pot_im(), 1000);
	
	cout <<"calculate simple_correlator on device..." << endl;
	usetimer noop;
	device.init_fermion_variables(&parameters, local_work_size, global_work_size, &inittimer);
	device.simple_correlator_device(&copytimer, &singletimer, &Mtimer, &scalarprodtimer, &latimer, &solvertimer, &dslashtimer, &Mdiagtimer,  local_work_size, global_work_size, 1000);
	device.finalize_fermions();
#endif

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
  
#else //_PERFORM_BENCHMARKS_

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////
	// Benchmarking
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////	
	
	benchmark_id = argv[2];

#ifndef _FERMIONS_

  int benchmarksteps1 = parameters.get_heatbathsteps();
	cout<<"perform HEATBATH-BENCHMARK with "<<benchmarksteps1<<" steps off each device operation..."<<endl;
	for(int i = 0; i<benchmarksteps1; i++){
		device.run_heatbath(parameters.get_beta(), local_work_size, global_work_size, &updatetime);
		device.run_overrelax(parameters.get_beta(), local_work_size, global_work_size, &overrelaxtime);
		device.gaugeobservables(local_work_size, global_work_size, &plaq, &tplaq, &splaq, &pol, &plaqtime, &polytime);
		device.get_gaugefield_from_device(gaugefield, &copytime);
	}
	
#else

	int benchmarksteps2 = parameters.get_thermalizationsteps();
	int cgmax = parameters.get_cgmax();
	cout<<"perform FERMION-BENCHMARK with "<<benchmarksteps2<<" steps off each device operation..."<<endl;
	
	//CP: set up testing field
	hmc_spinor_field in[SPINORFIELDSIZE];
	if(!use_eo){
		init_spinorfield_cold(in);
		device.copy_spinorfield_to_device(in, &copytimer);
	}
	else{
		//!!CP: this should be fine since only half the field is used but of course it is not nice...
		init_spinorfield_cold_eoprec(in);
		device.copy_eoprec_spinorfield_to_device(in, &copytimer);
	}
	device.init_fermion_variables(&parameters, local_work_size, global_work_size, &inittimer);
	device.perform_benchmark(benchmarksteps2, cgmax, local_work_size, global_work_size, &copytimer, &singletimer, &Mtimer, &scalarprodtimer, &latimer, &solvertimer, &dslashtimer, &Mdiagtimer);
	device.finalize_fermions();
	
#endif //_FERMIONS_
	
#endif //_PERFORM_BENCHMARKS_

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Final Output
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  totaltime.add();
#ifndef _PERFORM_BENCHMARKS_
  save_gaugefield(gaugefield, &parameters, nsteps);  
#endif
	time_output(
  	&totaltime, &inittime, &polytime, &plaqtime, &updatetime, &overrelaxtime, &copytime
#ifdef _FERMIONS_
, &inittimer, &singletimer, &Mtimer, &copytimer, &scalarprodtimer, &latimer, &solvertimer, &dslashtimer, &Mdiagtimer
#endif
	);
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // free variables
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  free(gaugefield);
  device.finalize();
  
  return HMC_SUCCESS;
}
