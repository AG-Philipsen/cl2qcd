#include "gaugefield_inversion.h"


hmc_error Gaugefield_inversion::init(int numdevs, cl_device_type* devicetypes, inputparameters* input_parameters, usetimer* timer){
	
	hmc_error err;
  
	//LZ: problem is finalize...
/*	hmc_gaugefield * gf_tmp = (hmc_gaugefield*) malloc(sizeof(hmc_gaugefield));
	err = set_gf (gf_tmp);*/
gf = (hmc_gaugefield*) malloc(sizeof(hmc_gaugefield));
	
	err = set_parameters (input_parameters);

	init_gaugefield(timer);

	err = set_num_ocl_devices (numdevs);

	if(numdevs != 1) {
		//LZ: so far, we only use !!! 1 !!! device
		//this needs generalisation to several devices and subsets!!!!!
		cerr << "only 1 device possible..." << endl;
	}

	if(get_num_ocl_devices() > 0){
		devices = new Opencl_fermions[get_num_ocl_devices()];
//		err = set_devices(devices_tmp);
	}


	for(int n = 0; n < get_num_ocl_devices(); n++) {
		cout << "init device #" << n << endl;
		(get_devices_f())[n].init(devicetypes[n], local_work_size, global_work_size, timer, get_parameters ());
		
	}

	return HMC_SUCCESS;
}

Opencl * Gaugefield_inversion::get_devices_f (){
  return  devices;
}
