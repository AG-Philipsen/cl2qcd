/** @file
 * Provides a class for gauge fields
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
#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include "exceptions.h"
#include "gaugefield_hybrid.h"
#include "opencl_module.h"
#include "opencl_module_ran.h"
#include "opencl_module_heatbath.h"

#include "logger.hpp"


/**
 * Version number.
 */
extern string const version;

/**
 * Class for the gaugefield. Includes initialization, device management for two devices doing a heatbath and tk calculation.
 *
 * @class Gaugefield
 */
class Gaugefield_heatbath : public Gaugefield_hybrid {
public:
	//init functions
	/**
	 * Initialize tasks
	 * Called by init()
	 */
	virtual void init_tasks();

	// proper finish
	/**
	 * Free variables
	 * Called by finalize()
	 */
	virtual void delete_variables();
	/**
	 * Free variables
	 * Called by finalize()
	 */
	virtual void finalize_opencl();

	// do the real work
	/**
	 * Performs 1 heatbath and nover overrelaxation steps.
	 * @param[in] nover number of overrelaxation steps per heatbath step
	 */
	void perform_tasks(int nover);

	// get methods
	/**
	 * Returns a pointer to heatbath task
	 */
	Opencl_Module_Heatbath* get_task_heatbath();

private:

	int task_heatbath;

};

#endif /* _GAUGEFIELDHEATBATHH_ */
