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

#include "exceptions.h"

#include "logger.hpp"


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
	 */
	void init(int numdevs, cl_device_type* devicetypes, inputparameters* input_parameters);

	/**
	 * Fills an array of device types according to inputparameters
	 */
	void init_devicetypes_array(cl_device_type* devicetypes, inputparameters* input_parameters);

	/**
	 * Initialize device types
	 *
	 * @todo This needs to be worked out in detail. So far it is assumed that numdevs[] has identical entries.
	 */
	void init(int* numdevs, int numdevtypes, cl_device_type* devicetypes, inputparameters* input_parameters);

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
	virtual void init_devices(cl_device_type* devicetypes);
	/**
	 * Initializes the gaugefield, to be called by init()
	 */
	void init_gaugefield();


	//communication
	/**
	 * Copy gaugefield to devices (here: to device).
	 */
	virtual void copy_gaugefield_to_devices();
	/**
	 * Copy gaugefield from devices (currently: from device) to host.
	 * @param[in,out] timer copy-time
	 */
	virtual void sync_gaugefield();
	/**
	 * Copy random array to devices (currently: to device).
	 */
	void copy_rndarray_to_devices();
	/**
	 * Copy random array from devices (currently: from device).
	 */
	void copy_rndarray_from_devices();

	//input/output, print, save functions!!
	/**
	 * Save gaugefield to file.
	 */
	void save(int number);
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
	 * Print gauge observables passed as function arguments, add iteration number to stdout.
	 * @param[in] plaq plaquette value
	 * @param[in] tplaq timelike plaquette value
	 * @param[in] splaq spatial plaquette value
	 * @param[in] pol Polyakov loop value
	 * @param[in] iter integer number that accompanies output
	 */
	void print_gaugeobservables(hmc_float plaq, hmc_float tplaq, hmc_float splaq, hmc_complex pol, int iter);

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
	 * @param[in] stdout print also to stdout
	 */
	virtual void print_gaugeobservables_from_devices(hmc_float * const plaq, hmc_float * const tplaq, hmc_float * const splaq, hmc_complex * const pol, const int i, const string gaugeoutname, int stdout);
	/**
	 * Print gauge observables calculated on device, add iteration number.
	 * @param[in,out] plaqtime time to calculate plaquette
	 * @param[in,out] polytime time to calculate Polyakov loop
	 * @param[in] i integer number that accompanies output
	 * @param[in] gaugeoutname name of output file
	 * @param[in] stdout print also to stdout
	 */
	virtual void print_gaugeobservables_from_devices(const int i, const string gaugeoutname, int stdout);


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
	Matrixsu3 * get_sgf ();

	/**
	 * Sets private member gaugefield u (structures)
	 */
	void set_sgf (Matrixsu3 * sgf_val);

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
	 */
	void set_devices (Opencl * devices_val);

	/**
	 * Sets private member * devices with index i
	 * @in index i
	 */
	void set_devices (Opencl * devices_val, int i);

	/**
	 * Returns private member num_ocl_devices
	 * @return num_ocl_devices
	 */
	int get_num_ocl_devices ();
	/**
	 * Sets private member num_ocl_devices
	 */
	void set_num_ocl_devices (int num);
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
	 */
	void set_num_device_types (int num);

	/**
	 * Sets private member * parameters
	 * @return parameters
	 */
	void set_parameters (inputparameters * parameters_val);

	/**
	 * Copies the gaugefield from pure array format to structure array format
	 * @param[in] gf pure array
	 * @param[out] sgf array of structs
	 */
	void copy_gaugefield_to_s_gaugefield (Matrixsu3 * sgfo, hmc_complex * gf);

	/**
	 * Copies the gaugefield from structure array format to pure array format
	 * @param[in] sgf array of structs
	 * @param[out] gf pure array
	 */
	void copy_s_gaugefield_to_gaugefield (hmc_complex * gf, Matrixsu3 * sgfo);

	/**
	 * Initializing the gaugefield consisting of structs for a hot start
	 * Not implemented yet, does a cold start
	 * @param field gaugefield
	 */
	void set_gaugefield_hot_new(Matrixsu3 * field);
	/**
	 * Initializing the gaugefield consisting of structs for a cold start
	 * @param field gaugefield
	 */
	void set_gaugefield_cold_new(Matrixsu3 * field);

	/**
	 * Set a value in the given gaugefield;
	 *
	 * @param[inout] field The field to set the value to.
	 * @param[in]    mu    Direction coordinate
	 * @param[in]    x     Space coordinate
	 * @param[in]    t     Time coordinate
	 * @param[in]    val   The value to set.
	 */
	void set_to_gaugefield(Matrixsu3 * field, const size_t mu, const size_t x, const size_t t, const Matrixsu3 val);

	/**
	 * Get a value from the given gaugefield;
	 *
	 * @param[in] field The field to set the value to.
	 * @param[in] mu    Direction coordinate
	 * @param[in] x     Space coordinate
	 * @param[in] t     Time coordinate
	 * @return The value stored at the given coordinates.
	 */
	Matrixsu3 get_from_gaugefield(const Matrixsu3 * field, const size_t mu, const size_t x, const size_t t) const;

	/**
	 * This method provides allocation for device double pointer
	 */
	void alloc_devicetypes();

	/**
	 * Returns private member * rndarray
	 * @return rndarray
	 */
	hmc_ocl_ran* get_rndarray();

	/*
	       * Returns numrndstates
	 * @return numrndstates
	 */
	size_t get_numrndstates();

	usetimer * get_copy_to();
	usetimer * get_copy_on();

protected:

	size_t get_num_hmc_gaugefield_elems();

private:


	Opencl ** devices;
	inputparameters* parameters;
	Matrixsu3 * sgf;

	int num_ocl_devices;
	int num_device_types;

	int numrndstates;
	size_t sizeof_rndarray;

	hmc_ocl_ran* rndarray;

	//bunch of timers
	//this is used to measure data-transfer to and from the device
	usetimer copy_to;
	//this is used to measure data-transfer on the device
	usetimer copy_on;
};

#endif /* _GAUGEFIELDH_ */
