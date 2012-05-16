/** @file
 * Provides a class for gauge fields
 *
 */
#ifndef _GAUGEFIELDHMCH_
#define _GAUGEFIELDHMCH_

#include <cstdlib>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "globaldefs.h"
#include "types.h"
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
	void md_update_gaugemomentum_fermion(hmc_float eps, usetimer * solvertimer, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void md_update_gaugemomentum_detratio(hmc_float eps, usetimer * solvertimer);
	void md_update_gaugefield(hmc_float eps);
	void init_gaugemomentum_spinorfield(usetimer * solvertimer);

	/**
	 * See hep-lat/0505020 for more infos on integrators.
	 * Currently, two Integrators are implemented: Leapfrog and 2MN.
	 * Leapfrog-Integration is the simplest second-order integration, thus the erros made are not
	 *  small with regard to 2MN-integration. However, it only requires one force-calculation per step.
	 * 2MN (2order Minimized Norm) needs two force-calculations per step, but has a much lower error then leapfrog. It also
	 *  features one tunable parameter, lambda.
	 * If one has the Hamiltonian:
	 *  H = p^2/2 + S(q) = T + V
	 * then one leapfrog step of length eps (= tau/number_of_steps) is
	 *  exp(eps(T + V) ) = exp(eps/2 T) exp( eps V ) exp( eps/2 T)
	 * and one 2MN-step is
	 *  exp(eps(T + V) ) = exp(lambda*eps T) exp( eps/2 V ) exp( (1 - 2lamdba) *eps T) exp( eps/2 V ) exp( lamdba*eps T)
	 * One can also use more than one timescale for the different pars to improve speed (e.g. hep-lat/0209037 and hep-lat/0506011v2).
	 * Up to now, 2 timescales are implemented (hardcoded). The Hamiltonian is split up
	 *  H = T + V_gauge + V_fermion = S_1 + S_2
	 * In case of the Leapfrog this means then
	 *  exp(eps(S_1) ) = exp(eps/2 T(V_gauge)) exp( eps V ) exp( eps/2 T (V_gauge) )
	 *  exp(eps(S_2) ) = exp(eps T(V_fermion))
	 *  -> exp(eps(T + V) ) = exp(eps/2 S_1) [ exp( eps/m S_2 ) ]^m exp( eps/2 S_1)
	 * In case of the 2MN this means then
	 *  exp(eps(S_1) ) = exp(lambda*eps T(V_gauge) ) exp( eps/2 V ) exp( (1 - 2lamdba) *eps T(V_gauge) ) exp( eps/2 V ) exp( lamdba*eps T(V_gauge) )
	 *  exp(eps(S_2) ) = exp(eps T(V_fermion))
	 *  -> exp(eps(T + V) ) = exp(lambda*eps/2 S_1) [ exp( eps/m S_2 ) ]^m exp( (1 - 2lamdba)*eps/2 S_1) [ exp( eps/m S_2 ) ]^m exp( lamdba*eps/2 S_1)
	 * In the program:
	 *  exp(eps T) == md_update_gaugemomentum(eps)
	 *  exp(eps V) == md_update_gaugefield(eps)
	 * NOTE:  Performing number_of_steps integrationsteps leads to "whole" steps for the momentum
	 *          ( exp(eps/2 T)exp(eps/2 T) = exp(eps T) )
	 */
	void integrator(usetimer * solvertimer);
	void leapfrog(usetimer * solvertimer);
	void twomn(usetimer * solvertimer);
	void calc_total_force(usetimer * solvertimer);
	void fermion_forces_call(usetimer * solvertimer, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void detratio_forces_call(usetimer * solvertimer);

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

#endif /* _GAUGEFIELDHMCH_ */
