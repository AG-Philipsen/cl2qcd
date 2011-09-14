/** @file
 * Basic OpenCL functionality
 */
#ifndef _OPENCLMODULECORRELATORH_
#define _OPENCLMODULECORRELATORH_

#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include "host_geometry.h"
#include "host_operations_complex.h"
#include "host_operations_gaugefield.h"
#include "globaldefs.h"
#include "types.h"
#include "host_use_timer.h"
#include "inputparameters.h"
#include "opencl_compiler.hpp"

#include "opencl_module.h"
#include "opencl_module_ran.h"
#include "opencl_module_spinors.h"

#include "exceptions.h"

/**
 * An OpenCL device
 *
 * This class wraps all operations on a device. Operations are always specific, e.g. each kernel and copy operation
 * has it's own wrapper function.
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Opencl_Module_Correlator : public Opencl_Module_Spinors {
public:

	// OpenCL specific methods needed for building/compiling the OpenCL program
	/**
	 * Collect the compiler options for OpenCL.
	 * Virtual method, allows to include more options in inherited classes.
	 */
	virtual void fill_collect_options(stringstream* collect_options);
	/**
	 * Collect the buffers to generate for OpenCL.
	 * Virtual method, allows to include more buffers in inherited classes.
	 */
	virtual void fill_buffers();
	/**
	 * Collect the kernels for OpenCL.
	 * Virtual method, allows to include more kernels in inherited classes.
	 */
	virtual void fill_kernels();
	/**
	 * Clear out the kernels,
	 * Virtual method, allows to clear additional kernels in inherited classes.
	 */
	virtual void clear_kernels();
	/**
	 * Clear out the buffers,
	 * Virtual method, allows to clear additional buffers in inherited classes.
	 */
	virtual void clear_buffers();

	/**
	 * comutes work-sizes for a kernel
	 * @todo autotune
	 * @param ls local-work-size
	 * @param gs global-work-size
	 * @param num_groups number of work groups
	 * @param dev_type type of device on which the kernel should be executed
	 * @param name name of the kernel for possible autotune-usage, not yet used!!
	 */
	virtual void get_work_sizes(const cl_kernel kernel, cl_device_type dev_type, size_t * ls, size_t * gs, cl_uint * num_groups);

	void create_point_source_device(cl_mem inout, int i, int spacepos, int timepos);

	void create_stochastic_source_device(cl_mem inout);

	void correlator_device(const cl_kernel correlator_kernel, cl_mem in, cl_mem correlator);

	/**
	 * Get kernel for correlator indicated by which
	 * @param[in] which string that identifies the correlator (ps or sc)
	 * @return correlator_kernel
	 */
	cl_kernel get_correlator_kernel(string which);

	/////////////////////////////////////////////////
	//functions to get private variables
	cl_mem get_clmem_corr();
	cl_mem get_clmem_source();

#ifdef _PROFILING_
	//CP: if PROFILING is activated, one needs a timer for each kernel

	usetimer timer_create_point_source;
	usetimer timer_create_stochastic_source;

	usetimer timer_ps_correlator;

	/**
	 * Return the timer connected to a specific kernel.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual usetimer* get_timer(char * in);

	/**
	 * Return amount of bytes read and written by a specific kernel per call.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual int get_read_write_size(char * in, inputparameters * parameters);

	/**
	 * Print the profiling information to a file.
	 *
	 * @param filename Name of file where data is appended.
	 * @param parameters inputparameters
	 */
	void virtual print_profiling(std::string filename);
#endif

private:
	////////////////////////////////////
	//kernels

	cl_kernel create_point_source;
	cl_kernel create_stochastic_source;

	//Observables
	cl_kernel correlator_ps;
	cl_kernel correlator_sc;

	cl_mem clmem_source;

	cl_mem clmem_corr;

protected:
	ClSourcePackage basic_correlator_code;

};

#endif //OPENCLMODULECORRELATORH
