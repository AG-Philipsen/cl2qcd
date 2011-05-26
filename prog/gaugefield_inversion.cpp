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

hmc_error Gaugefield_inversion::free_devices(){
  if(get_num_ocl_devices() > 0)
    delete [] get_devices_fermions();
  return HMC_SUCCESS;
}

Opencl_fermions * Gaugefield_inversion::get_devices_fermions (){
  return  (Opencl_fermions*)get_devices();
}
