#include "gaugefield_inversion.h"


void Gaugefield_inversion::init_devices(cl_device_type* devicetypes)
{
	if(get_num_ocl_devices() != 1) {
		//LZ: so far, we only use !!! 1 !!! device
		//this needs generalisation to several devices and subsets!!!!!
		cerr << "only 1 device possible..." << endl;
	}

	if(get_num_ocl_devices() > 0) {
		Opencl_fermions* dev_tmp = new Opencl_fermions[get_num_ocl_devices()];
		alloc_devicetypes();
		set_devices(dev_tmp);
	}

	for(int n = 0; n < get_num_ocl_devices(); n++) {
		logger.trace() << "init device #" << n ;
		get_devices_fermions()[n].init(devicetypes[n], get_parameters(),this->get_numrndstates());
	}
	return;
}

void Gaugefield_inversion::finalize()
{
	Gaugefield::finalize();
	for(int n = 0; n < get_num_ocl_devices(); n++)
		get_devices_fermions()[n].finalize_fermions();
	return;
}

void Gaugefield_inversion::free_devices()
{
	if(get_num_ocl_devices() > 0)
		delete [] get_devices_fermions();
	return;
}

Opencl_fermions * Gaugefield_inversion::get_devices_fermions ()
{
	return  (Opencl_fermions*)get_devices();
}

void Gaugefield_inversion::perform_inversion_pointsource_ps_corr_devices(usetimer* solvertimer){
	int use_eo = get_parameters()->get_use_eo();

  if(use_eo==false){
    get_devices_fermions()[0].set_spinorfield_cold_device(get_devices_fermions()[0].get_clmem_inout());
  }
  else{
    get_devices_fermions()[0].set_eoprec_spinorfield_cold_device(get_devices_fermions()[0].get_clmem_inout_eoprec());		
  }

  get_devices_fermions()[0].set_zero_spinorfield_device(get_devices_fermions()[0].get_clmem_corr());

  for(int k=0; k<12; k++) {
    if(use_eo == false){
      get_devices_fermions()[0].create_point_source_device(get_devices_fermions()[0].get_clmem_source(), k,0,0);
      get_devices_fermions()[0].solver_device(M_call, get_devices_fermions()[0].get_clmem_inout(), get_devices_fermions()[0].get_clmem_source(), get_devices_fermions()[0].get_clmem_gaugefield(), solvertimer, get_parameters()->get_cgmax());
      //CP: add solution to former ones...
			get_devices_fermions()[0].saxpy_device(get_devices_fermions()[0].get_clmem_inout(), get_devices_fermions()[0].get_clmem_corr(), get_devices_fermions()[0].get_clmem_minusone(), get_devices_fermions()[0].get_clmem_corr());
    }
    else{
      get_devices_fermions()[0].create_point_source_eoprec_device(get_devices_fermions()[0].get_clmem_source_even(), get_devices_fermions()[0].get_clmem_source_odd(), get_devices_fermions()[0].get_clmem_gaugefield(), k,0,0);
      get_devices_fermions()[0].solver_eoprec_device(Aee_call, get_devices_fermions()[0].get_clmem_inout(), get_devices_fermions()[0].get_clmem_inout_eoprec(), get_devices_fermions()[0].get_clmem_source_even(), get_devices_fermions()[0].get_clmem_source_odd(), get_devices_fermions()[0].get_clmem_gaugefield(), solvertimer, get_parameters()->get_cgmax());
      //CP: add solution to former ones... This is the same call as without eoprec since the eoprec solver saves the normal field also in clmem_inout!!
      get_devices_fermions()[0].saxpy_device(get_devices_fermions()[0].get_clmem_inout(), get_devices_fermions()[0].get_clmem_corr(), get_devices_fermions()[0].get_clmem_minusone(), get_devices_fermions()[0].get_clmem_corr());
    }
  }
  return;
}
