#include "gaugefield_hmc.h"


hmc_error Gaugefield_hmc::init_devices(cl_device_type* devicetypes, usetimer* timer)
{
	if(get_num_ocl_devices() != 1) {
		//LZ: so far, we only use !!! 1 !!! device
		//this needs generalisation to several devices and subsets!!!!!
		cerr << "only 1 device possible..." << endl;
	}

	if(get_num_ocl_devices() > 0) {
		Opencl_hmc* dev_tmp = new Opencl_hmc[get_num_ocl_devices()];
		set_devices(dev_tmp);
	}


	for(int n = 0; n < get_num_ocl_devices(); n++) {
		cout << "init device #" << n << endl;
		get_devices_hmc()[n].init(devicetypes[n], timer, get_parameters());
	}
	return HMC_SUCCESS;
}

hmc_error Gaugefield_hmc::finalize()
{
	hmc_error err = HMC_SUCCESS;
	err |= Gaugefield_hmc::finalize();
	for(int n = 0; n < get_num_ocl_devices(); n++)
		err |= get_devices_hmc()[n].finalize_hmc();
	return err;
}

hmc_error Gaugefield_hmc::free_devices()
{
	if(get_num_ocl_devices() > 0)
		delete [] get_devices_hmc();
	return HMC_SUCCESS;
}

Opencl_hmc * Gaugefield_hmc::get_devices_hmc ()
{
	return  (Opencl_hmc*)get_devices();
}

/*
hmc_error Gaugefield_hmc::perform_inversion_pointsource_ps_corr_devices(usetimer* copytimer, usetimer* singletimer, usetimer* Mtimer, usetimer* scalarprodtimer, usetimer* latimer, usetimer* dslashtimer, usetimer* Mdiagtimer, usetimer* solvertimer){

  Gaugefield_inversion.perform_inversion_pointsource_ps_corr_device(copytimer, singletimer, Mtimer, scalarprodtimer, latimer, dslashtimer, Mdiagtimer, solvertimer);

  return HMC_SUCCESS;
}
*/
