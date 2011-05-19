/** @file
 * Provides a class for gauge fields
 *
 * @todo: Apply the class also to the hmc, benchmark and inverter
 *
 */
#ifndef _GAUGEFIELDH_
#define _GAUGEFIELDH_

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
#include "host_input.h"
#include "host_readgauge.h"
#include "host_writegaugefield.h"
#include "host_use_timer.h"
#include "opencl.h"
#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif


/**
 * Version number.
 */
extern string const version;

/**
 * Class for the gaugefield. Includes initialization, device management and heatbath update.
 *
 * @class Gaugefield
 * @todo: Also incorporate hmc and solver capabilities. Generalize to several devices.
 */
class Gaugefield {
public:
	/**
	 * Initialize gaugefield and devices.
	 *
	 * @param[in] numdevs Number of wanted devices (so far, only 1 makes sense).
	 * @param[in] devicetypes Array of wanted cl_device_types for the devices.
	 * @param[in] input_parameters instance of inputparameters that contains information from input file
	 * @param[in,out] timer Return initialization time.
	 * @return Error code as defined in hmcerrs.h
	 */
        virtual hmc_error init(int numdevs, cl_device_type* devicetypes, inputparameters* input_parameters, usetimer* timer);
	/**
	 * Free gaugefield and device allocations.
	 */
	hmc_error finalize();
	/**
	 * Save gaugefield to file.
	 */
	hmc_error save(int number);
	/**
	 * Copy gaugefield to devices (currently: to device).
	 * @param[in,out] timer copy-time
	 */
	hmc_error copy_gaugefield_to_devices(usetimer* timer);
	/**
	 * Copy random array to devices (currently: to device).
	 * @param[in,out] timer copy-time
	 */
	hmc_error copy_rndarray_to_devices(hmc_rndarray host_rndarray,  usetimer* timer);
	/**
	 * Copy random array from devices (currently: from device).
	 * @param[in,out] timer copy-time
	 */
	hmc_error copy_rndarray_from_devices(hmc_rndarray rndarray, usetimer* timer);
	/**
	 * Copy gaugefield from devices (currently: from device) to host.
	 * @param[in,out] timer copy-time
	 */
	hmc_error sync_gaugefield(usetimer* timer);
	/**
	 * Print gauge observables calculated from host gaugefield to stdout.
	 * @param[in,out] timer time to calculate plaquette
	 * @param[in,out] timer2 time to calculate Polyakov loop
	 */
	void print_gaugeobservables(usetimer* timer, usetimer * timer2);
	/**
	 * Print gauge observables calculated from host gaugefield to stdout, add iteration number.
	 * @param[in,out] timer time to calculate plaquette
	 * @param[in,out] timer2 time to calculate Polyakov loop
	 * @param[in] iter Integer number that accompanies output
	 */
	void print_gaugeobservables(usetimer* timer, usetimer * timer2, int iter);
	/**
	 * Print gauge observables calculated from host gaugefield to file, add iteration number.
	 * @param[in,out] timer time to calculate plaquette
	 * @param[in,out] timer2 time to calculate Polyakov loop
	 * @param[in] iter integer number that accompanies output
	 * @param[in] file name of output file
	 */
	void print_gaugeobservables(usetimer* timer, usetimer * timer2, int iter, std::string file);
	/**
	 * Print gauge observables passed as function arguments, add iteration number.
	 * @param[in] plaq plaquette value
	 * @param[in] tplaq timelike plaquette value
	 * @param[in] splaq spatial plaquette value
	 * @param[in] pol Polyakov loop value
	 * @param[in] iter integer number that accompanies output
	 * @param[in] file name of output file
	 */
	void print_gaugeobservables(hmc_float plaq, hmc_float tplaq, hmc_float splaq, hmc_complex pol, int iter, std::string filename);
	/**
	 * Print gauge observables calculated on device, add iteration number, return values to program.
	 * @param[in,out] plaq pointer to plaquette value
	 * @param[in,out] tplaq pointer to timelike plaquette value
	 * @param[in,out] splaq pointer to spatial plaquette value
	 * @param[in,out] pol pointer to Polyakov loop value
	 * @param[in,out] plaqtime time to calculate plaquette
	 * @param[in,out] polytime time to calculate Polyakov loop
	 * @param[in] i integer number that accompanies output
	 * @param[in] gaugeoutname name of output file
	 * @return Error code as defined in hmcerrs.h
	 */
	hmc_error print_gaugeobservables_from_devices(hmc_float * const plaq, hmc_float * const tplaq, hmc_float * const splaq, hmc_complex * const pol, usetimer * const plaqtime, usetimer * const polytime, const int i, const string gaugeoutname);
	/**
	 * Print gauge observables calculated on device, add iteration number.
	 * @param[in,out] plaqtime time to calculate plaquette
	 * @param[in,out] polytime time to calculate Polyakov loop
	 * @param[in] i integer number that accompanies output
	 * @param[in] gaugeoutname name of output file
	 * @return Error code as defined in hmcerrs.h
	 */
	hmc_error print_gaugeobservables_from_devices(usetimer * const plaqtime, usetimer * const polytime, const int i, const string gaugeoutname);

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
	/**
	 * Returns private member gaugefield
	 * @return The gaugefield
	 */
	hmc_gaugefield * get_gf ();
	
	/**
	 * Sets private member gaugefield
	 * @return Error code as defined in hmcerrs.h
	 */
	hmc_error set_gf (hmc_gaugefield * gf_val);

	/**
	 * Returns private member * devices
	 * @return devices
	 */
	Opencl*  get_devices ();
	
	/**
	 * Sets private member * devices
	 * @return Error code as defined in hmcerrs.h
	 */
	hmc_error set_devices (Opencl * devices_val);
	
	/**
	 * Returns private member num_ocl_devices
	 * @return num_ocl_devices
	 */
	int get_num_ocl_devices ();
	
	/**
	 * Sets private member num_ocl_devices
	 * @return Error code as defined in hmcerrs.h
	 */
	hmc_error set_num_ocl_devices (int num);
	
	/**
	 * Returns private member * parameters
	 * @return parameters
	 */
	inputparameters * get_parameters ();
	
	/**
	 * Sets private member * parameters
	 * @return parameters
	 */
	hmc_error set_parameters (inputparameters * parameters_val);
	
	/**
	 * ?
	 * @param[in,out] params ?
	 */
	void print_info_source(sourcefileparameters* params);
	
	/**
	 * Initializes the gaugefield
	 * @param[in,out] timer timer for initialization
	 * @return Error code as defined in hmcerrs.h
	 */
	hmc_error init_gaugefield(usetimer* timer);
	
	//gaugeobservables, on host!!
	hmc_float plaquette(hmc_float* tplaq, hmc_float* splaq);
	hmc_float plaquette();
	hmc_complex polyakov();
	hmc_complex spatial_polyakov(int dir);
	
	hmc_gaugefield * gf;	
	Opencl * devices;
	//Die sollten auch von den Folgeklassen aus zugreifbar sein
	inputparameters* parameters;
	
private:
	int num_ocl_devices;
// 	inputparameters* parameters;

};

#endif /* _GAUGEFIELDH_ */
