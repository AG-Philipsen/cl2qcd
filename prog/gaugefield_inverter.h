/** @file
 * Provides a class for gauge fields
 *
 */
#ifndef _GAUGEFIELDINVERTERH_
#define _GAUGEFIELDINVERTERH_

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
#include "opencl_module.h"
#include "opencl_module_fermions.h"
#include "opencl_module_correlator.h"
#include "gaugefield_hybrid.h"

#include "logger.hpp"


/**
 * Version number.
 */
extern string const version;

/**
 * Class for the gaugefield. Includes initialization, device management for multiple devices. Performs calculation of correlators.
 *
 * @class Gaugefield_inverter
 */
class Gaugefield_inverter : public Gaugefield_hybrid {
 public:

	/**
	 * Initialize class.
	 * Helper function called by init()
	 * Virtual to allow proper modification in inherited classes
	 */
	virtual void init_tasks();

	/**
	 * Free gaugefield and device allocations.
	 * Called by finalize().
	 * Virtual to allow proper modification in inherited classes
	 */
	virtual void delete_variables();
	/**
	 * Free gaugefield and device allocations.
	 * Called by finalize().
	 * Virtual to allow proper modification in inherited classes
	 */
	virtual void finalize_opencl();


	// the real job:
	/**
	 * Perform the inversion and store result to solution_buffer
	 * @param[in,out] solver_timer
	 */
	void perform_inversion(usetimer* solver_timer);

	/**
	 * Calculate flavour multiplet correlators from private solution_buffer and store them to a file
	 * @param[in] corr_fn filename
	 */
	void flavour_doublet_correlators(string corr_fn);

	// get methods
	/**
	 * Returns a pointer to solver task
	 */
	Opencl_Module_Fermions* get_task_solver();

	/**
	 * Returns a pointer to correlator task
	 */
	Opencl_Module_Correlator* get_task_correlator();

protected:

private:

	spinorfield* solution_buffer;

	int task_solver;
	int task_correlator;

};

#endif /* _GAUGEFIELDINVERTERH_ */
