/** @file
 * Provides a class for gauge fields
 *
 */
#ifndef _GAUGEFIELDHEATBATHKAPPAH_
#define _GAUGEFIELDHEATBATHKAPPAH_

#include <cstdlib>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "globaldefs.h"
#include "types.h"
#include "host_operations_gaugefield.h"
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
#include "opencl_module_kappa.h"

#include "logger.hpp"


/**
 * Class for the gaugefield. Includes initialization, device management for two devices doing a heatbath and tk calculation.
 *
 * @class Gaugefield
 */
class Gaugefield_heatbath_kappa : public Gaugefield_hybrid {
public:
	/**
	 * Create a new kappa heatbath gaugefield.
	 *
	 * \param params The input parameters of the application
	 */
	Gaugefield_heatbath_kappa(const meta::Inputparameters& params)
		: Gaugefield_hybrid(params) { };

	//init functions
	/**
	 * Initialize tasks
	 * Called by init()
	 */
	virtual void init_tasks() override;

	// proper finish
	/**
	 * Free variables
	 * Called by finalize()
	 */
	virtual void delete_variables() override;
	/**
	 * Free variables
	 * Called by finalize()
	 */
	virtual void finalize_opencl() override;

	// do the real work
	/**
	 * Perform the heatbath only (e.g. for thermalization)
	 * @param[in] nheat number of heatbath steps
	 * @param[in] nover number of overrelaxation steps per heatbath step
	 */
	void perform_heatbath(int nheat, int nover);
	/**
	 * Perform all the tasks
	 * @param[in] nheat number of heatbath steps
	 * @param[in] nover number of overrelaxation steps per heatbath step
	 */
	void perform_tasks(int nheat, int nover);
	/**
	 * Perform all the tasks, perform tuning
	 * @param[in] nheat number of heatbath steps
	 * @param[in] nover number of overrelaxation steps per heatbath step
	 * @param[out] nheat_optimal optimal number of heatbath steps for further iterations
	 */
	void perform_tasks(int nheat, int nover, int* nheat_optimal);
	/**
	 * Print kappa coefficient to screen, add iteration number
	 * @param[in] iter iteration number
	 */
	void print_kappa(int iter);
	/**
	 * Print kappa coefficient to file, add iteration number
	 * @param[in] iter iteration number
	 * @param[in] filename name of file to which the output line is to be appended
	 */
	void print_kappa(int iter, std::string filename);

	// get methods
	/**
	 * Returns a pointer to heatbath task
	 */
	Opencl_Module_Heatbath* get_task_heatbath();
	/**
	 * Returns a pointer to kappa task
	 */
	Opencl_Module_Kappa* get_task_kappa();

private:

	int task_heatbath;
	int task_kappa;


};

#endif /* _GAUGEFIELDHEATBATHKAPPAH_ */
