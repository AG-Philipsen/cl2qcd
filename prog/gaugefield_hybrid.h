/** @file
 * Provides a class for gauge fields
 *
 */
#ifndef _GAUGEFIELDHYBRIDH_
#define _GAUGEFIELDHYBRIDH_

#include <cstdlib>
#include <fstream>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "hardware/system.hpp"

#include "globaldefs.h"
#include "types.h"
#include "host_operations_gaugefield.h"
#include "host_readgauge.h"
#include "host_writegaugefield.h"
#include "host_use_timer.h"

#include "exceptions.h"
#include "opencl_module.h"

#include "logger.hpp"


/**
 * Version number.
 *
 * @deprecated move this into some specific header or so
 */
extern std::string const version;

/**
 * Class for the gaugefield. Includes initialization, device management for multiple devices.
 *
 * @class Gaugefield
 */
class Gaugefield_hybrid {
public:

	//constructor is empty
	/**
	 * Create a new gaugefield.
	 *
	 * @param[in] params points to an instance of inputparameters
	 */
	Gaugefield_hybrid(const hardware::System * system)
		: system(system), parameters(system->get_inputparameters()) { };

	//init functions
	/**
	 * Initialize class.
	 * @param[in] numtasks sets the number of different tasks
	 * @param[in] primary_device_type defines which device type (CPU/GPU) is to be preferred
	 *
	 * @deprecated move to constructor
	 */
	void init(int numtasks, cl_device_type primary_device_type);
	/**
	 * Initialize class.
	 * Helper function called by init()
	 */
	void init_opencl();
	/**
	 * Initialize devices.
	 * Helper function called by init_opencl()
	 * @param[in] ndev number of device to be initialized
	 */
	void init_devices(int ndev);

	/**
	 * Initialize class.
	 * Helper function called by init()
	 * Virtual to allow proper modification in inherited classes
	 */
	virtual void init_tasks();

	// proper finish
	/**
	 * Free gaugefield and device allocations.
	 * Called by destructor.
	 */
	void finalize();
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


	// gaugefield specific methods
	/**
	 * Initializes the gaugefield, called by init()
	 * calls init_gaugefield(sourcefile) with the name from inputparameters
	 */
	void init_gaugefield();
	/**
	 * Initializes the gaugefield, called by init()
	 * @param name of gaugefield file
	 */
	void init_gaugefield(const char*);

	/**
	 * Initializing the gaugefield consisting of structs for a hot start
	 * Not implemented yet, throws an exception
	 * @param field gaugefield
	 */
	void set_gaugefield_hot(Matrixsu3 * field);
	/**
	 * Initializing the gaugefield consisting of structs for a cold start
	 * @param field gaugefield
	 */
	void set_gaugefield_cold(Matrixsu3 * field);

	/**
	 * Set the number of devices available
	 * @param num number of devices
	 */
	void set_num_devices(int num);
	/**
	 * Get the number of devices available
	 * @return num_devices
	 */
	int get_num_devices();


	// communication

	/**
	 * Distribute the host copy of gaugefield to all devices
	 */
	void copy_gaugefield_to_all_tasks();
	/**
	 * Distribute the host copy of gaugefield to a specific devices
	 * @param[in] ntask target task
	 */
	void copy_gaugefield_to_task(int ntask);
	/**
	 * Copy the gaugefield from task with number ntask to host
	 * @param[in] ntask target task
	 */
	void copy_gaugefield_from_task(int ntask);
	/**
	 * OpenCL barrier.
	 * Copy gaugefield from ntask_reference to host and distribute to other tasks
	 * @param[in] ntask target task
	 */
	void synchronize(int ntask_reference);

	// get and set methods

	/**
	 * Sets private member num_tasks
	 * @param[in] num number to be set
	 */
	void set_num_tasks (int num);
	/**
	 * Returns private member num_tasks
	 * @return num_tasks
	 */
	int get_num_tasks ();
	/**
	 * Returns private member * parameters
	 * @return parameters
	 */
	const meta::Inputparameters& get_parameters ();
	/**
	 * Sets private member gaugefield u (structures)
	 */
	void set_sgf (Matrixsu3 * sgf_val);
	/**
	 * Returns pointer to gaugefield u (structures)
	 * @return The gaugefield
	 */
	Matrixsu3 * get_sgf ();
	/**
	 * Get device assigned to given task
	 * @param[in] ntask id for task
	 * @return device cl_device_id for given task
	 */
	hardware::Device* get_device_for_task(int ntask);

	// output methods
	/**
	 * Create a name for a gaugefield configuration file
	 * @param[in] number number of the file
	 */
	std::string create_configuration_name(int number);
	/**
	 * Save gaugefield to a file with name conf.number
	 * @param[in] number number to be added to file name
	 */
	void save(int number);
	/**
	 * Save gaugefield to a file with given name
	 * @param[in] outputfile name of file
	 */
	void save(std::string outputfile);
	/**
	 * Return plaquette value (calculated from host gaugefield)
	 */
	hmc_float plaquette();
	/**
	 * Return plaquette value (calculated from host gaugefield)
	 * @param[out] tplaq timelike plaquette
	 * @param[out] splaq spatial plaquette
	 */
	hmc_float plaquette(hmc_float* tplaq, hmc_float* splaq);
	/**
	 * Return value of Polyakov loop (calculated from host gaugefield)
	 */
	hmc_complex polyakov();
	/**
	 * Return value of spatial Polyakov loop (calculated from host gaugefield)
	 * @param[in] dir spatial direction for Polyakov loop calculation
	 */
	hmc_complex spatial_polyakov(int dir);
	/**
	 * Print gaugeobservables to screen, add iteration number
	 * @param[in] iter iteration number
	 */
	void print_gaugeobservables(int iter);
	/**
	 * Print gaugeobservables to file, add iteration number
	 * @param[in] iter iteration number
	 * @param[in] filename name of file to which the output line is to be appended
	 */
	void print_gaugeobservables(int iter, std::string filename);
	/**
	 * Print gaugeobservables from task to screen, add iteration number
	 * @param[in] iter iteration number
	 * @param[in] ntask id of target task
	 */
	void print_gaugeobservables_from_task(int iter, int ntask);
	/**
	 * Print gaugeobservables from task to file, add iteration number
	 * @param[in] iter iteration number
	 * @param[in] ntask id of target task
	 * @param[in] filename name of file to which the output line is to be appended
	 */
	void print_gaugeobservables_from_task(int iter, int ntask, std::string filename);

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
	 * Get the number of elements in the gaugefield.
	 */
	size_t get_num_gaugefield_elems() const;

	/**
	* Create the gaugefield from an array of floats as used by by ILDG.
	*
	* @param[out] gaugefield Pointer to the new storage location.
	* @param[in] gaugefield_tmp Field in IDLG format
	* @param[in] check Size of the ILDG field.
	* @todo Replace hmc_gaugefield type by s_gaugefield type (LZ)
	*/
	void copy_gaugefield_from_ildg_format(Matrixsu3 * gaugefield, hmc_float * gaugefield_tmp, int check);
	/**
	* Create the IDLG representation of the given gaugefield.
	*
	* @param[out] dest The location to store the ILDG representation to
	* @param[in] source The gaugefield in the internal representation
	* @todo Replace hmc_gaugefield type by s_gaugefield type (LZ)
	*/
	void copy_gaugefield_to_ildg_format(hmc_float * dest, Matrixsu3 * source);

	void print_profiling(std::string filename);

	/**
	 * Returns parameters read from sourcefile
	 */
	const sourcefileparameters get_parameters_source();

	/**
	 * This checks the parameters from a sourcefile
	 * against the inputparameters.
	 * It indicates if minor parameters mismatch
	 * and throws an exception if major parameters do
	 */
	void check_sourcefileparameters();

protected:
	cl_device_type* devicetypes;
	Opencl_Module ** opencl_modules;

private:
	/**
	 * Initialize class.
	 * Helper function called by init()
	 */
	void init_devicetypearray(cl_device_type primary_device_type);

	const hardware::System * system;
	const meta::Inputparameters& parameters;
	sourcefileparameters parameters_source;
	Matrixsu3 * sgf;
	int num_tasks;
	int num_devices;
	int* device_id_for_task;

	//OpenCL:
	cl_device_id* cl_devices;
	hardware::Device** devices;
};

#endif /* _GAUGEFIELDHYBRIDH_ */
