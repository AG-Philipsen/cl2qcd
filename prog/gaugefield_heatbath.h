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
#include "hmcerrs.h"
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
	virtual hmc_error finalize();
	/**
	 * Free device, called by finalize
	 */
	virtual hmc_error free_devices();
	/**
	 * Initializes the devices, to be called by init()
	 * @return Error code as defined in hmcerrs.h
	 * @param devicetypes array of cl_device_type handles
	 * @param[in,out] timer timer for initialization
	 */
	virtual	hmc_error init_devices(cl_device_type* devicetypes);


	//calculations on device
	/**
	 * Perform a number of heatbath and (afterwards) overrelaxation steps.
	 * @param[in] nheat number of heatbath steps
	 * @param[in] nover number of overrelaxation steps
	 * @param[in,out] timer_heat time for heatbath steps
	 * @param[in,out] timer_over time for overrelaxation steps
	 */
	hmc_error heatbath(const int nheat, const int nover, usetimer * const timer_heat, usetimer * const timer_over);
	/**
	 * Perform a number of heatbath steps.
	 * @param[in] nheat number of heatbath steps
	 * @param[in,out] timer_heat time for heatbath steps
	 */
	hmc_error heatbath(const int nheat, usetimer * const timer_heat);
	/**
	 * Perform one heatbath step.
	 * @param[in,out] timer time for heatbath step
	 */
	hmc_error heatbath(usetimer * const timer);
	/**
	 * Perform one overrelaxation step.
	 * @param[in,out] timer time for overrelaxation step
	 */
	hmc_error overrelax(usetimer * const timer);

	//access to private members
	
	/**
	 * Returns private member * devices
	 * @return devices
	 */
	Opencl_heatbath*  get_devices_heatbath ();

private:

};

#endif /* _GAUGEFIELDHHEATBATH_ */
