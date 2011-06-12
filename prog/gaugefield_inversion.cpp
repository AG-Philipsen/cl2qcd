#include "gaugefield_inversion.h"


hmc_error Gaugefield_inversion::init_devices(cl_device_type* devicetypes, usetimer* timer){
  if(get_num_ocl_devices() != 1) {
    //LZ: so far, we only use !!! 1 !!! device
    //this needs generalisation to several devices and subsets!!!!!
    cerr << "only 1 device possible..." << endl;
  }
  
  if(get_num_ocl_devices() > 0) {
    Opencl_fermions* dev_tmp = new Opencl_fermions[get_num_ocl_devices()];
    set_devices(dev_tmp);
  }

  
  for(int n = 0; n < get_num_ocl_devices(); n++) {
	  cout << "init device #" << n << endl;
	  get_devices_fermions()[n].init(devicetypes[n], timer, get_parameters());	  
  }
  return HMC_SUCCESS;
}

hmc_error Gaugefield_inversion::finalize(){
  hmc_error err = HMC_SUCCESS;
  err |= Gaugefield::finalize();
  for(int n = 0; n < get_num_ocl_devices(); n++)
    err |= get_devices_fermions()[n].finalize_fermions();	  
  return err;
}

hmc_error Gaugefield_inversion::free_devices(){
  if(get_num_ocl_devices() > 0)
    delete [] get_devices_fermions();
  return HMC_SUCCESS;
}

Opencl_fermions * Gaugefield_inversion::get_devices_fermions (){
  return  (Opencl_fermions*)get_devices();
}

//CP: this is deprecated 
/*
hmc_error Gaugefield_inversion::perform_inversion_pointsource_ps_corr_host(){
  //CP: one needs bicgstab here for M
  int use_cg = FALSE;
	
  int use_eo = get_parameters()->get_use_eo();

  hmc_spinor_field in[SPINORFIELDSIZE];
  hmc_spinor_field phi[SPINORFIELDSIZE];
  init_spinorfield_cold(in);

  //pseudo scalar, flavour multiplet
  hmc_complex correlator_ps[NSPACE];
  for(int z=0; z<NSPACE; z++) {
    correlator_ps[z].re = 0;
    correlator_ps[z].im = 0;
  }
  cout << "start iterations of solver.. " << endl;
 for(int k=0; k<NC*NSPIN; k++) {
    if(use_eo == FALSE){
      hmc_spinor_field b[SPINORFIELDSIZE];
      cout << "create pointsource.." << endl;
      create_point_source(get_parameters(),k,0,0,b);
      cout << "start solver.." << endl;
      solver(get_parameters(), in, b, get_sgf(), use_cg, phi);
    cout<<"solver done" << endl;
    }
    else{
      hmc_eoprec_spinor_field be[EOPREC_SPINORFIELDSIZE];
      hmc_eoprec_spinor_field bo[EOPREC_SPINORFIELDSIZE];
      
      create_point_source_eoprec(get_parameters(), k,0,0, get_gf(), be,bo);
      solver_eoprec(get_parameters(), in, be, bo, get_sgf(), use_cg, phi);
      
    }
    for(int timepos = 0; timepos<NTIME; timepos++) {
      for(int spacepos = 0; spacepos<VOLSPACE; spacepos++) {
	for(int alpha = 0; alpha<NSPIN; alpha++) {
	  for(int c = 0; c<NC; c++) {
	    // int j = spinor_element(alpha,c);
	    int n = spinor_field_element(alpha, c, spacepos, timepos);
	    int z = get_spacecoord(spacepos, 3);
	    hmc_complex tmp = phi[n];
	    hmc_complex ctmp = complexconj(&tmp);
	    hmc_complex incr = complexmult(&ctmp,&tmp);
	    correlator_ps[z].re += incr.re;
	    correlator_ps[z].im += incr.im;
	  }
	}
      }
    }
  }
  
  printf("pseudo scalar correlator:\n");
  for(int z=0; z<NSPACE; z++) {
    printf("%d\t(%e,%e)\n",z,correlator_ps[z].re,correlator_ps[z].im);
  }

  return HMC_SUCCESS;
}
*/

hmc_error Gaugefield_inversion::perform_inversion_pointsource_ps_corr_devices(usetimer* copytimer, usetimer* singletimer, usetimer* Mtimer, usetimer* scalarprodtimer, usetimer* latimer, usetimer* dslashtimer, usetimer* Mdiagtimer, usetimer* solvertimer){

  //this uses a BiCGStab inverter on device


  //global and local work sizes;
  //LZ: should eventually be moved inside opencl_fermions class
#ifdef _USE_GPU_
	const size_t ls = NUM_THREADS; /// @todo have local work size depend on kernel properties (and device? autotune?)
#else
	const size_t ls = 1; // nothing else makes sens on CPU
#endif

#ifdef _USE_GPU_
	size_t gs = 4 * NUM_THREADS * get_devices()[0].max_compute_units; /// @todo autotune
#else
	size_t gs = get_devices()[0].max_compute_units;
#endif

	const cl_uint num_groups = (gs + ls - 1) / ls;
	gs = ls * num_groups;

  /** @todo here one has to introduce more timer instead of noop*/
  usetimer noop;

  int use_eo = get_parameters()->get_use_eo();

  cout << "calc simple_propagator on the device..." << endl;
  
  if(use_eo==FALSE){
    get_devices_fermions()[0].set_spinorfield_cold_device(ls, gs, &noop);
    //get_devices_fermions()[0].copy_spinorfield_to_device(in, copytimer);
  }
  else{
    //!!CP: this should be fine since only half the field is used but of course it is not nice...
    //init_spinorfield_cold_eoprec(in);
    //get_devices_fermions()[0].copy_eoprec_spinorfield_to_device(in, copytimer);		
  }

  for(int k=0; k<NC*NSPIN; k++) {
    if(use_eo == FALSE){
      cout << "create point source" << endl;
      get_devices_fermions()[0].create_point_source_device(k,0,0,ls, gs, latimer);
      cout << "start solver " << endl;
      get_devices_fermions()[0].solver_device(copytimer, singletimer, Mtimer, scalarprodtimer, latimer, dslashtimer, Mdiagtimer, solvertimer, ls, gs, get_parameters()->get_cgmax());
    }
    else{
      //get_devices_fermions()[0].create_point_source_eoprec_device(k,0,0,ls, gs, latimer, dslashtimer, Mdiagtimer);
      //get_devices_fermions()[0].solver_eoprec_device(phi, copytimer, singletimer, Mtimer, scalarprodtimer, latimer, dslashtimer, Mdiagtimer, solvertimer, ls, gs, get_parameters()->get_cgmax());
    }
  }

  get_devices_fermions()[0].ps_correlator_device(ls, gs, &noop);

  //copy spinorfield back to host (this should not be done anymore in the end)
  //get_devices_fermions()[0].get_spinorfield_from_device(phi, copytimer);
  
  return HMC_SUCCESS;
}
