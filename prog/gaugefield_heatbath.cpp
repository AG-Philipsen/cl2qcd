#include "gaugefield_heatbath.h"

#include "logger.hpp"

hmc_error Gaugefield_heatbath::init_devices(cl_device_type* devicetypes)
{
// 	if(get_num_ocl_devices() != 1) {
// 		//LZ: so far, we only use !!! 1 !!! device
// 		//this needs generalisation to several devices and subsets!!!!!
// 		logger.error() << "only 1 device possible...";
// 	}

	if(get_num_ocl_devices() > 0) {
		Opencl_heatbath* dev_tmp = new Opencl_heatbath[get_num_ocl_devices()];
		alloc_devicetypes();
		set_devices(dev_tmp);
	}


	for(int n = 0; n < get_num_ocl_devices(); n++) {
		logger.debug() << "init device #" << n;
		get_devices_heatbath()[n].init(devicetypes[n], get_parameters());
	}
	return HMC_SUCCESS;
}


hmc_error Gaugefield_heatbath::heatbath(usetimer * const timer)
{
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...
	hmc_error err = get_devices_heatbath()[0].run_heatbath(get_parameters()->get_beta(), timer);
	return err;
}

hmc_error Gaugefield_heatbath::overrelax(usetimer * const timer)
{
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...

	hmc_error err = get_devices_heatbath()[0].run_overrelax(get_parameters()->get_beta(), timer);
	return err;
}

hmc_error Gaugefield_heatbath::heatbath(const int nheat, const int nover, usetimer * const timer_heat, usetimer * const timer_over)
{
	hmc_error err = HMC_SUCCESS;
	for(int i = 0; i < nheat; i++) err |= heatbath(timer_heat);
	for(int i = 0; i < nover; i++) err |= overrelax(timer_over);
	return err;
}

hmc_error Gaugefield_heatbath::heatbath(const int nheat, usetimer * const timer_heat)
{
	hmc_error err = HMC_SUCCESS;
	for(int i = 0; i < nheat; i++) err |= heatbath(timer_heat);
	return err;
}

hmc_error Gaugefield_heatbath::free_devices()
{
	if(get_num_ocl_devices() > 0)
		delete [] get_devices_heatbath();
	return HMC_SUCCESS;
}

Opencl_heatbath * Gaugefield_heatbath::get_devices_heatbath ()
{
  return  (Opencl_heatbath*)get_devices();
}

hmc_error Gaugefield_heatbath::finalize()
{
	hmc_error err = HMC_SUCCESS;
	err |= Gaugefield::finalize();
	return HMC_SUCCESS;
}
