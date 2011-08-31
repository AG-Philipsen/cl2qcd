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

#include "logger.hpp"


/**
 * Version number.
 */
extern string const version;

/**
 * Class for the gaugefield. Includes initialization, device management for multiple devices.
 *
 * @class Gaugefield
 */
class Gaugefield_hybrid {
public:
  //init/finalize functions

  /**
   * Initialize class.
   */
  void init(int numtasks, cl_device_type primary_device_type, inputparameters* input_parameters);

  void init_devicetypearray(cl_device_type primary_device_type);
  void init_opencl();
  void init_devices();
  void init_random_arrays();

  /**
   * Free gaugefield and device allocations.
   */
  void finalize();

  void delete_variables();
  void finalize_opencl();

  /**
   * Initializes the gaugefield, to be called by init()
   */
  void init_gaugefield();
  
  /**
   * Initializing the gaugefield consisting of structs for a hot start
   * Not implemented yet, does a cold start
   * @param field gaugefield
   */
  void set_gaugefield_hot(s_gaugefield * field);

  /**
   * Initializing the gaugefield consisting of structs for a cold start
   * @param field gaugefield
   */
  void set_gaugefield_cold(s_gaugefield * field);

  void copy_gaugefield_to_all_devices();
  void copy_gaugefield_to_device(int ntask);
  void copy_gaugefield_from_device(int ntask);
  void synchronize(int ntask_reference);

  
  /**
   * Sets private member num_tasks
   */
  void set_num_tasks (int num);

  /**
   * Returns private member num_tasks
   * @return num_tasks
   */
  int get_num_tasks ();

  /**
   * Sets private member * parameters
   * @return parameters
   */
  void set_parameters (inputparameters * parameters_val);

  /**
   * Returns private member * parameters
   * @return parameters
   */
  inputparameters * get_parameters ();

  /**
   * Sets private member gaugefield u (structures)
   */
  void set_sgf (s_gaugefield * sgf_val);

  /**
   * Returns pointer to gaugefield u (structures)
   * @return The gaugefield
   */
  s_gaugefield * get_sgf ();

  /**
   * Returns private member * rndarray for given task
   * @param[in] ntask task identifier
   * @param[in] ndevice device identifier
   * @return rndarray
   */
  hmc_ocl_ran* get_rndarray(int ntask);

  /**
   * Returns device type for given task.
   */
  cl_device_type get_device_type(int ntask);

  /**
   * Copies the gaugefield from pure array format to structure array format
   * @param[in] gf pure array
   * @param[out] sgf array of structs
   */
  void copy_gaugefield_to_s_gaugefield (s_gaugefield * sgfo, hmc_gaugefield * gf);
  
  /**
   * Copies the gaugefield from structure array format to pure array format
   * @param[in] sgf array of structs
   * @param[out] gf pure array
   */
  void copy_s_gaugefield_to_gaugefield (hmc_gaugefield * gf, s_gaugefield * sgfo);

  void save(int number);
  void save(string outputfile);

  hmc_float plaquette();
  hmc_float plaquette(hmc_float* tplaq, hmc_float* splaq);
  hmc_complex polyakov();
  hmc_complex spatial_polyakov(int dir);

  void print_gaugeobservables(int iter);
  void print_gaugeobservables(int iter, std::string filename);

private:

  inputparameters* parameters;
  s_gaugefield * sgf;
  int num_tasks;
  
  int* numrndstates;
  size_t* sizeof_rndarray;
  hmc_ocl_ran** rndarray;
	
  cl_device_type* devicetypes;
  //	Opencl ** devices;
  
  //OpenCL:
  cl_platform_id platform;
  cl_context context;
  cl_command_queue* queue;
  cl_device_id* devices;

  string* device_double_extension;
  cl_uint* max_compute_units;

  cl_mem clmem_gaugefield;

};

#endif /* _GAUGEFIELDHYBRIDH_ */
