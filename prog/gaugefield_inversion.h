/** @file
 * Provides class Gaugefield_inversion, inherited from Gaugefield, adds solver capabilities
 *
 */
#ifndef _GAUGEFIELDINVERSIONH_
#define _GAUGEFIELDINVERSIONH_

#include "gaugefield.h"
#include "opencl_fermions.h"

class Gaugefield_inversion : public Gaugefield {
  public:
	 /**
	 * Initialize gaugefield and devices for fermion matrix inversion
	 *
	 * @param[in] numdevs Number of wanted devices (so far, only 1 makes sense).
	 * @param[in] devicetypes Array of wanted cl_device_types for the devices.
	 * @param[in] input_parameters instance of inputparameters that contains information from input file
	 * @param[in,out] timer Return initialization time.
	 * @return Error code as defined in hmcerrs.h
	 */
	hmc_error init(int numdevs, cl_device_type* devicetypes, inputparameters* input_parameters, usetimer* timer);
	
	/**
	 * Returns private member Opencl_fermions * devices
	 * @return devices of type Opencl_fermions
	 */
	Opencl * get_devices_f();
	
	//	Opencl_fermions * devices;
};


#endif //_GAUGEFIELDINVERSIONH_
