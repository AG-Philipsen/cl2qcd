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
  //init/finalize functions
	/**
	 * Initialize gaugefield and devices.
	 *
	 * @param[in] numdevs Number of wanted devices (so far, only 1 makes sense).
	 * @param[in] devicetypes Array of wanted cl_device_types for the devices.
	 * @param[in] input_parameters instance of inputparameters that contains information from input file
	 * @param[in,out] timer Return initialization time.
	 * @return Error code as defined in hmcerrs.h
	 */
        hmc_error init(int numdevs, cl_device_type* devicetypes, inputparameters* input_parameters, usetimer* timer);

	/**
	 * Initialize gaugefield and different device types
	 *
	 * @todo This needs to be worked out in detail. So far it is assumed that numdevs[] has identical entries.
	 */
        hmc_error init(int* numdevs, int numdevtypes, cl_device_type* devicetypes, inputparameters* input_parameters, usetimer* timer);
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
	virtual	hmc_error init_devices(cl_device_type* devicetypes, usetimer* timer);
	/**
	 * Initializes the gaugefield, to be called by init()
	 * @param[in,out] timer timer for initialization
	 * @return Error code as defined in hmcerrs.h
	 */
	hmc_error init_gaugefield(usetimer*  timer);
	

	//communication
	/**
	 * Copy gaugefield to devices (here: to device).
	 * @param[in,out] timer copy-time
	 */
	virtual hmc_error copy_gaugefield_to_devices(usetimer* timer);
	/**
	 * Copy gaugefield from devices (currently: from device) to host.
	 * @param[in,out] timer copy-time
	 */
	virtual hmc_error sync_gaugefield(usetimer* timer);


	//input/output, print, save functions!!
	/**
	 * Save gaugefield to file.
	 */
	hmc_error save(int number);
	/**
	 * Print information from source file (if any)
	 * @param[in,out] params instance of sourcefileparameters
	 */
	void print_info_source(sourcefileparameters* params);
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
	virtual hmc_error print_gaugeobservables_from_devices(hmc_float * const plaq, hmc_float * const tplaq, hmc_float * const splaq, hmc_complex * const pol, usetimer * const plaqtime, usetimer * const polytime, const int i, const string gaugeoutname);
	/**
	 * Print gauge observables calculated on device, add iteration number.
	 * @param[in,out] plaqtime time to calculate plaquette
	 * @param[in,out] polytime time to calculate Polyakov loop
	 * @param[in] i integer number that accompanies output
	 * @param[in] gaugeoutname name of output file
	 * @return Error code as defined in hmcerrs.h
	 */
	virtual hmc_error print_gaugeobservables_from_devices(usetimer * const plaqtime, usetimer * const polytime, const int i, const string gaugeoutname);
	

	//gaugeobservables, on host!!
	/**
	 * Calculate plaquette on host.
	 * @param[out] tplaq timelike plaquette
	 * @param[out] splaq spatial plaquette
	 * @return plaquette
	 * @todo we still need to transform to the old hmc-format, we should implement new operations for host code
	 */
	hmc_float plaquette(hmc_float* tplaq, hmc_float* splaq);
	/**
	 * Calculate plaquette on host.
	 * @return plaquette
	 */
	hmc_float plaquette();
	/**
	 * Calculate Polyakov loop (in time direction) on host.
	 * @return Polyakov loop
	 */
	hmc_complex polyakov();
	/**
	 * Calculate Polyakov loop (in spatial direction dir) on host.
	 * @param[in] dir spatial direction to be used, dir=1,2,3
	 * @return Polyakov loop
	 */
	hmc_complex spatial_polyakov(int dir);


	//access to private members
	/**
	 * Returns pointer to gaugefield u (structures)
	 * @return The gaugefield
	 */
	s_gaugefield * get_sgf ();
	
	/**
	 * Sets private member gaugefield u (structures)
	 * @return Error code as defined in hmcerrs.h
	 */
	hmc_error set_sgf (s_gaugefield * sgf_val);
	
	/**
	 * Returns private member * devices
	 * @return devices
	 */
	Opencl*  get_devices ();

	/**
	 * Returns private member * devices with index i
	 * @in index of device type
	 * @return devices
	 */
	Opencl*  get_devices (int i);

	/**
	 * Sets private member * devices
	 * @return Error code as defined in hmcerrs.h
	 */
	hmc_error set_devices (Opencl * devices_val);

	/**
	 * Sets private member * devices with index i
	 * @in index i
	 * @return Error code as defined in hmcerrs.h
	 */
	hmc_error set_devices (Opencl * devices_val, int i);

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
	 * Returns private member num_device_types
	 * @return num_device_types
	 */
	int get_num_device_types ();
	/**
	 * Sets private member num_device_types
	 * @return Error code as defined in hmcerrs.h
	 */
	hmc_error set_num_device_types (int num);

	/**
	 * Sets private member * parameters
	 * @return parameters
	 */
	hmc_error set_parameters (inputparameters * parameters_val);
	
	/**
	 * Copies the gaugefield from pure array format to structure array format
	 * @param[in] gf pure array
	 * @param[out] sgf array of structs
	 * @return Error code as defined in hmcerrs.h
	 */
	hmc_error copy_gaugefield_to_s_gaugefield (s_gaugefield * sgfo, hmc_gaugefield * gf);
	
	/**
	 * Copies the gaugefield from structure array format to pure array format
	 * @param[in] sgf array of structs
	 * @param[out] gf pure array
	 * @return Error code as defined in hmcerrs.h
	 */
	hmc_error copy_s_gaugefield_to_gaugefield (hmc_gaugefield * gf, s_gaugefield * sgfo);
	
	/**
	 * Initializing the gaugefield consisting of structs for a hot start
	 * Not implemented yet, does a cold start
	 * @param field gaugefield
	 * @return Error code as defined in hmcerrs.h
	 */
	hmc_error set_gaugefield_hot_new(s_gaugefield * field);
	/**
	 * Initializing the gaugefield consisting of structs for a cold start
	 * @param field gaugefield
	 * @return Error code as defined in hmcerrs.h
	 */
	hmc_error set_gaugefield_cold_new(s_gaugefield * field);
	
	/**
	 * This method provides allocation for device double pointer
	 * @return Error code as defined in hmcerrs.h
	 */
	hmc_error alloc_devicetypes();

private:
	Opencl ** devices;
	inputparameters* parameters;
	s_gaugefield * sgf;
	
	int num_ocl_devices;
	int num_device_types;
};

#endif /* _GAUGEFIELDH_ */
