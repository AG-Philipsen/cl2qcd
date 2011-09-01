#include "gaugefield_heatbath.h"

#include "logger.hpp"

void Gaugefield_heatbath::init_devices(cl_device_type* devicetypes)
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
		get_devices_heatbath()[n].init(devicetypes[n], get_parameters(),this->Gaugefield::get_numrndstates());
	}
	return;
}


void Gaugefield_heatbath::heatbath()
{
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...
	get_devices_heatbath()[0].run_heatbath();
	return;
}

void Gaugefield_heatbath::overrelax()
{
	//LZ: so far, we only use !!! 1 !!! device
	// this function needs to be generalised to several devices and definition of subsets...

	get_devices_heatbath()[0].run_overrelax();
	return;
}

void Gaugefield_heatbath::heatbath(const int nheat, const int nover)
{
	for(int i = 0; i < nheat; i++) heatbath();
	for(int i = 0; i < nover; i++) overrelax();
	return;
}

void Gaugefield_heatbath::heatbath(const int nheat)
{
	for(int i = 0; i < nheat; i++) heatbath();
	return;
}

void Gaugefield_heatbath::free_devices()
{
	if(get_num_ocl_devices() > 0)
		delete [] get_devices_heatbath();
	return;
}

Opencl_heatbath * Gaugefield_heatbath::get_devices_heatbath ()
{
  return  (Opencl_heatbath*)get_devices();
}

void Gaugefield_heatbath::finalize()
{
	Gaugefield::finalize();
	return;
}
