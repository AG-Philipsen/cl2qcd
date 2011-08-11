#include "gaugefield_inversion.h"


hmc_error Gaugefield_inversion::init_devices(cl_device_type* devicetypes)
{
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
		get_devices_fermions()[n].init(devicetypes[n], get_parameters());
	}
	return HMC_SUCCESS;
}

hmc_error Gaugefield_inversion::finalize()
{
	hmc_error err = HMC_SUCCESS;
	err |= Gaugefield::finalize();
	for(int n = 0; n < get_num_ocl_devices(); n++)
		err |= get_devices_fermions()[n].finalize_fermions();
	return err;
}

hmc_error Gaugefield_inversion::free_devices()
{
	if(get_num_ocl_devices() > 0)
		delete [] get_devices_fermions();
	return HMC_SUCCESS;
}

Opencl_fermions * Gaugefield_inversion::get_devices_fermions ()
{
	return  (Opencl_fermions*)get_devices();
}

hmc_error Gaugefield_inversion::perform_inversion_pointsource_ps_corr_devices(usetimer* copytimer, usetimer* singletimer, usetimer* Mtimer, usetimer* scalarprodtimer, usetimer* latimer, usetimer* dslashtimer, usetimer* Mdiagtimer, usetimer* solvertimer){
	//this uses a BiCGStab inverter on device

	//global and local work sizes;
	//LZ: should eventually be moved inside opencl_fermions class
#ifdef _USEGPU_
	const size_t ls = NUMTHREADS; /// @todo have local work size depend on kernel properties (and device? autotune?)
#else
	const size_t ls = 1; // nothing else makes sens on CPU
#endif

#ifdef _USEGPU_
	size_t gs = 4 * NUMTHREADS * get_devices()[0].max_compute_units; /// @todo autotune
#else
	size_t gs = get_devices()[0].max_compute_units;
#endif

	const cl_uint num_groups = (gs + ls - 1) / ls;
	gs = ls * num_groups;

  /** @todo here one has to introduce more timer instead of noop*/
  usetimer noop;

	int use_eo = get_parameters()->get_use_eo();

  if(use_eo==FALSE){
    get_devices_fermions()[0].set_spinorfield_cold_device(ls, gs, &noop);
  }
  else{
    get_devices_fermions()[0].set_eoprec_spinorfield_cold_device(ls, gs, &noop);		
  }

  get_devices_fermions()[0].set_correlator_field_zero_device(ls, gs, latimer);

  for(int k=0; k<12; k++) {
    if(use_eo == FALSE){
      get_devices_fermions()[0].create_point_source_device(k,0,0,ls, gs, latimer);
      get_devices_fermions()[0].solver_device(copytimer, singletimer, Mtimer, scalarprodtimer, latimer, dslashtimer, Mdiagtimer, solvertimer, ls, gs, get_parameters()->get_cgmax());
      //CP: add solution to former ones...
      get_devices_fermions()[0].add_solution_to_correlator_field_device(ls, gs, latimer);
    }
    else{
      get_devices_fermions()[0].create_point_source_eoprec_device(k,0,0,ls, gs, latimer, dslashtimer, Mdiagtimer);
      get_devices_fermions()[0].solver_eoprec_device(copytimer, singletimer, Mtimer, scalarprodtimer, latimer, dslashtimer, Mdiagtimer, solvertimer, ls, gs, get_parameters()->get_cgmax());
      //CP: add solution to former ones... This is the same call as without eoprec since the eoprec solver saves the normal field also in clmem_inout!!
      get_devices_fermions()[0].add_solution_to_correlator_field_device(ls, gs, latimer);
    }
  }
	
	//CP: this should be called from outside...
  //get_devices_fermions()[0].ps_correlator_device(ls, gs, &noop);

  return HMC_SUCCESS;
}
