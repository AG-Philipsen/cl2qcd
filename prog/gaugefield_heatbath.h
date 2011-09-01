/** @file
 * Provides a class for gauge fields including heatbath operations
 *
 */
#ifndef _GAUGEFIELDHEATBATHH_
#define _GAUGEFIELDHEATBATHH_

#include <cstdlib>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "globaldefs.h"
#include "types.h"
#include "host_operations_complex.h"
#include "host_operations_gaugefield.h"
#include "inputparameters.h"
#include "host_readgauge.h"
#include "host_writegaugefield.h"
#include "host_use_timer.h"
#include "opencl.h"
#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include "gaugefield.h"
#include "opencl_heatbath.h"

/**
 * Class for the gaugefield. Adds heatbath update to the standard Gaugefield class.
 *
 * @class Gaugefield_heatbath
 */
class Gaugefield_heatbath : public Gaugefield {
public:
  //init/finalize functions

	/**
	 * Free gaugefield and device allocations.
	 */
	virtual void finalize();
	/**
	 * Free device, called by finalize
	 */
	virtual void free_devices();
	/**
	 * Initializes the devices, to be called by init()
	 * @param devicetypes array of cl_device_type handles
	 * @param[in,out] timer timer for initialization
	 */
	virtual	void init_devices(cl_device_type* devicetypes);


	//calculations on device
	/**
	 * Perform a number of heatbath and (afterwards) overrelaxation steps.
	 * @param[in] nheat number of heatbath steps
	 * @param[in] nover number of overrelaxation steps
	 */
	void heatbath(const int nheat, const int nover);
	/**
	 * Perform a number of heatbath steps.
	 * @param[in] nheat number of heatbath steps
	 */
	void heatbath(const int nheat);
	/**
	 * Perform one heatbath step.
	 */
	void heatbath();
	/**
	 * Perform one overrelaxation step.
	 */
	void overrelax();

	//access to private members
	
	/**
	 * Returns private member * devices
	 * @return devices
	 */
	Opencl_heatbath*  get_devices_heatbath ();

private:

};

#endif /* _GAUGEFIELDHHEATBATH_ */
