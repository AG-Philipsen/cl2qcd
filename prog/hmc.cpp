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

#ifdef _FERMIONS_

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
	cout<<"perform FERMION-BENCHMARK with "<<benchmarksteps2<<" steps off each device operation..."<<endl;
	
	//CP: set up testing field
	hmc_spinor_field in[SPINORFIELDSIZE];
	if(!use_eo){
		init_spinorfield_cold(in);
	}
	else{
		//!!CP: this should be fine since only half the field is used but of course it is not nice...
		init_spinorfield_cold_eoprec(in);
	}
	hmc_float resid;
	device.init_fermion_variables(&parameters, local_work_size, global_work_size, &inittimer);
	for(int i = 0; i<benchmarksteps2; i++){
		if(!use_eo){
			copy_spinorfield_to_device(in, copytimer);
			convert_to_kappa_format_device(clmem_inout, ls, gs, latimer);
			set_zero_spinorfield_device(clmem_v, localsize, globalsize, &latimer);
			saxpy_device(clmem_rn, clmem_source, clmem_one, clmem_rn, localsize, globalsize, latimer);
			copy_complex_device(clmem_one, clmem_alpha, singletimer);
			set_complex_to_scalar_product_device(clmem_rhat, clmem_rn, clmem_rho_next, local_work_size, global_work_size, scalarprodtimer);
			set_complex_to_ratio_device(clmem_rho_next, clmem_rho, clmem_tmp1, singletimer);
			saxsbypz_device(clmem_p, clmem_v, clmem_rn, clmem_beta, clmem_tmp2, clmem_p, local_work_size, global_work_size, latimer);
			M_device(clmem_inout, clmem_rn, localsize, globalsize, Mtimer, dslashtimer, Mdiagtimer);
			set_float_to_global_squarenorm_device(clmem_rn, clmem_resid, local_work_size, global_work_size, scalarprodtimer);
			copy_float_from_device(clmem_resid, &resid, copytimer);
			copy_spinor_device(clmem_rn, clmem_rhat, singletimer);
			set_complex_to_product_device(clmem_tmp1, clmem_tmp2, clmem_beta, singletimer);
			create_point_source_device(k,0,0,ls, gs, latimer);
			solver_device(phi, copytimer, singletimer, Mtimer, scalarprodtimer, latimer, dslashtimer, Mdiagtimer, solvertimer, ls, gs, cgmax);
			convert_from_kappa_format_device(clmem_inout, clmem_inout, ls, gs, latimer);
			get_spinorfield_from_device(out, copytimer);
		}
		else{
			copy_eoprec_spinorfield_to_device(in, copytimer);
			set_zero_spinorfield_eoprec_device(clmem_v_eoprec, localsize, globalsize, latimer); 
			Aee_device(clmem_inout_eoprec, clmem_rn_eoprec, localsize, globalsize, Mtimer, singletimer, dslashtimer, Mdiagtimer, latimer);
			copy_eoprec_spinor_device(clmem_rn_eoprec, clmem_rhat_eoprec, singletimer);
			copy_complex_device(clmem_one, clmem_alpha, singletimer);
			set_complex_to_ratio_device(clmem_alpha, clmem_omega, clmem_tmp2, singletimer);
			set_complex_to_product_device(clmem_tmp1, clmem_tmp2, clmem_beta, singletimer);
			set_complex_to_scalar_product_eoprec_device(clmem_rhat_eoprec, clmem_rn_eoprec, clmem_rho_next, local_work_size, global_work_size, scalarprodtimer);
			saxsbypz_eoprec_device(clmem_p_eoprec, clmem_v_eoprec, clmem_rn_eoprec, clmem_beta, clmem_tmp2, clmem_p_eoprec, local_work_size, global_work_size, latimer);
			saxpy_eoprec_device(clmem_t_eoprec, clmem_s_eoprec, clmem_omega, clmem_rn_eoprec, local_work_size, global_work_size, latimer);
			set_float_to_global_squarenorm_eoprec_device(clmem_rn_eoprec, clmem_resid, local_work_size, global_work_size, scalarprodtimer);
			copy_float_from_device(clmem_resid, &resid, copytimer);
			create_point_source_eoprec_device(k,0,0,ls, gs, latimer, dslashtimer, Mdiagtimer);
			solver_eoprec_device(phi, copytimer, singletimer, Mtimer, scalarprodtimer, latimer, dslashtimer, Mdiagtimer, solvertimer, ls, gs, cgmax);
			convert_to_kappa_format_eoprec_device(clmem_inout_eoprec, ls, gs, latimer);
			convert_from_kappa_format_eoprec_device(clmem_inout_eoprec, clmem_inout_eoprec, ls, gs, latimer);
			get_eoprec_spinorfield_from_device(phi_even, copytimer);
		}
	}
	device.finalize_fermions();
	
#endif //_FERMIONS_
	
#endif //_PERFORM_BENCHMARKS_

  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  // Final Output
  //////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  totaltime.add();
  save_gaugefield(gaugefield, &parameters, benchmarksteps1);  
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
