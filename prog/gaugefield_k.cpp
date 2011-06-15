#include "gaugefield_k.h"

hmc_error Gaugefield_k::init_devices(cl_device_type* devicetypes, usetimer* timer){
	if(get_num_ocl_devices() != 1) {
		//LZ: so far, we only use !!! 1 !!! device
		//this needs generalisation to several devices and subsets!!!!!
		cerr << "only 1 device possible..." << endl;
	}

	if(get_num_ocl_devices() > 0) {
	  Opencl_k* dev_tmp = new Opencl_k[get_num_ocl_devices()];
	  set_devices(dev_tmp);
	}


	for(int n = 0; n < get_num_ocl_devices(); n++) {
		cout << "init device #" << n << endl;
		get_devices_k()[n].init(devicetypes[n], timer, get_parameters());
		
	}
	return HMC_SUCCESS;
}

hmc_error Gaugefield_k::free_devices(){
  if(get_num_ocl_devices() > 0)
    delete [] get_devices_k();
  return HMC_SUCCESS;
}

hmc_float Gaugefield_k::get_kappa_karsch (){
	return kappa_karsch_val;
}
	
hmc_float Gaugefield_k::get_kappa_clover (){
	return kappa_clover_val;
}

hmc_error Gaugefield_k::set_kappa_karsch (hmc_float in){
	kappa_karsch_val = in;
	return HMC_SUCCESS;
}

hmc_error Gaugefield_k::set_kappa_clover (hmc_float in){
	kappa_clover_val = in;
	return HMC_SUCCESS;
  
}

hmc_error Gaugefield_k::kappa_karsch_gpu (usetimer* timer){
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...
	hmc_error err;
 	err = get_devices_k()[0].run_kappa_karsch_gpu(get_parameters()->get_beta(), timer, &kappa_karsch_val);
	return err;
}

hmc_error Gaugefield_k::kappa_clover_gpu (usetimer* timer){
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...
	hmc_error err;
 	err = get_devices_k()[0].run_kappa_clover_gpu(get_parameters()->get_beta(), timer, &kappa_clover_val);
	return err;
}

Opencl_k * Gaugefield_k::get_devices_k (){
  return  (Opencl_k*)get_devices();
}
