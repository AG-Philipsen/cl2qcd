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
#include "opencl_module.h"
#include "opencl_module_fermions.h"
#include "opencl_module_correlator.h"
#include "gaugefield_hybrid.h"

#include "logger.hpp"

/**
 * Class for the gaugefield. Includes initialization, device management for multiple devices. Performs calculation of correlators.
 *
 * @class Gaugefield_inverter
 */
class Gaugefield_inverter : public Gaugefield_hybrid {
public:
	/**
	 * Create a new inverter gaugefield.
	 *
	 * \param params The input parameters of the application
	 */
	Gaugefield_inverter(const hardware::System * system)
		: Gaugefield_hybrid(system),
		  clmem_corr(nullptr),
		  clmem_source_solver(nullptr),
		  clmem_source_corr(nullptr)
	{ };

	/**
	 * Initialize class.
	 * Helper function called by init()
	 * Virtual to allow proper modification in inherited classes
	 */
	virtual void init_tasks() override;

	/**
	 * Free gaugefield and device allocations.
	 * Called by finalize().
	 * Virtual to allow proper modification in inherited classes
	 */
	virtual void delete_variables() override;
	/**
	 * Free gaugefield and device allocations.
	 * Called by finalize().
	 * Virtual to allow proper modification in inherited classes
	 */
	virtual void finalize_opencl() override;

	/**
	 * Invert fermionmatrix M corresponding to the upper flavour in Nf=2 LQCD.
	 * @param[in, out] inout spinorfield with initial guess and solution
	 * @param[in] source source to invert on
	 * @param[in] gf gaugefield configuration to invert on
	 * @param[in] solvertimer measures time needed for the inversion
	 */
	void invert_M_nf2_upperflavour(const hardware::buffers::Plain<spinor> * inout, const hardware::buffers::Plain<spinor> * source, const hardware::buffers::SU3 * gf, usetimer * solvertimer);

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
	void flavour_doublet_correlators(std::string corr_fn);

	// get methods
	/**
	 * Returns a pointer to solver task
	 */
	Opencl_Module_Fermions* get_task_solver();

	/**
	 * Returns a pointer to correlator task
	 */
	Opencl_Module_Correlator* get_task_correlator();

	void create_sources();

	const hardware::buffers::Plain<spinor> * get_clmem_corr();
	const hardware::buffers::Plain<spinor> * get_clmem_source_corr();

protected:

private:

	const hardware::buffers::Plain<spinor> * clmem_corr;
	const hardware::buffers::Plain<spinor> * clmem_source_solver;
	const hardware::buffers::Plain<spinor> * clmem_source_corr;

	spinor* solution_buffer;
	spinor* source_buffer;

	int task_solver;
	int task_correlator;

};

#endif /* _GAUGEFIELDINVERTERH_ */
