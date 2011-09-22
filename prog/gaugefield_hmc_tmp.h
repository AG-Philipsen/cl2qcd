/** @file
 * Provides a class for gauge fields
 *
 */
#ifndef _GAUGEFIELDHMCTMPH_
#define _GAUGEFIELDHMCTMPH_

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
#include "opencl_module_hmc.h"
#include "gaugefield_hybrid.h"
#include "types_hmc.h"

#include "logger.hpp"


/**
 * Version number.
 */
extern string const version;

/**
 * Class for the gaugefield. Includes initialization, device management for multiple devices.
 *
 * @class Gaugefield_inverter
 */
class Gaugefield_hmc : public Gaugefield_hybrid {
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
	void perform_hmc_step(hmc_observables *obs, int iter, hmc_float rnd_number, usetimer* solver_timer);
	
	void print_hmcobservables(hmc_observables obs, int iter, std::string filename);
	void print_hmcobservables(hmc_observables obs, int iter);
	
	void md_update_gaugemomentum(hmc_float eps, usetimer * solvertimer);
	void md_update_gaugemomentum_gauge(hmc_float eps);
	void md_update_gaugemomentum_fermion(hmc_float eps, usetimer * solvertimer);
	void md_update_gaugefield(hmc_float eps);
	void init_gaugemomentum_spinorfield();
	
	void integrator(usetimer * solvertimer);
	void leapfrog(usetimer * solvertimer);
	void twomn(usetimer * solvertimer);
	void calc_total_force(usetimer * solvertimer);
	void calc_gauge_force();
	void calc_fermion_force(usetimer * solvertimer);
	
	// get methods
	/**
	 * Returns a pointer to solver task
	 * @param[in] dev Number of the device one wants to choose
	 */
	Opencl_Module_Hmc* get_task_hmc(int dev);
	
protected:

private:

	int task_hmc;

};

#endif /* _GAUGEFIELDHMCTMPH_ */
