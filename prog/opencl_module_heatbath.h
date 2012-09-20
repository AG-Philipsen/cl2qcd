/** @file
 * Heatbath for OpenCL
 */
#ifndef _OPENCLMODULEHEATBATHH_
#define _OPENCLMODULEHEATBATHH_

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
#include "host_operations_gaugefield.h"
#include "globaldefs.h"
#include "types.h"
#include "host_use_timer.h"
#include "host_random.h"
#include "opencl_compiler.hpp"

#include "opencl_module.h"
#include "opencl_module_ran.h"

#include "exceptions.h"

/**
 * An OpenCL device
 *
 * Adds heatbath to basic Opencl_Module class
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Opencl_Module_Heatbath : public Opencl_Module_Ran {
public:
	/**
	 * Empty constructor.
	 *
	 * @param[in] params points to an instance of inputparameters
	 */
	Opencl_Module_Heatbath(const meta::Inputparameters& params, hardware::Device * device)
		: Opencl_Module_Ran(params, device) {};
	/**
	 * Collect the compiler options for OpenCL.
	 * Virtual method, allows to include more options in inherited classes.
	 */
	virtual void fill_collect_options(std::stringstream* collect_options) override;

	/**
	 * Collect the buffers to generate for OpenCL.
	 * Virtual method, allows to include more buffers in inherited classes.
	 */
	virtual void fill_buffers() override;

	/**
	 * Collect the kernels for OpenCL.
	 * Virtual method, allows to include more kernels in inherited classes.
	 */
	virtual void fill_kernels() override;

	/**
	 * Clear out the kernels,
	 * Virtual method, allows to clear additional kernels in inherited classes.
	 */
	virtual void clear_kernels() override;

	/**
	 * Clear out the buffers,
	 * Virtual method, allows to clear additional buffers in inherited classes.
	 */
	virtual void clear_buffers() override;

	/**
	 * Perform one heatbath step.
	 */
	void run_heatbath();

	/**
	 * Perform one overrelaxation step.
	 */
	void run_overrelax();

	/**
	 * Add specific work_size determination for this child class
	 */
	virtual void get_work_sizes(const cl_kernel kernel, cl_device_type dev_type, size_t * ls, size_t * gs, cl_uint * num_groups) override;

protected:

private:

	cl_kernel heatbath_odd;
	cl_kernel heatbath_even;
	cl_kernel overrelax_odd;
	cl_kernel overrelax_even;


#ifdef _PROFILING_
	//CP: if PROFILING is activated, one needs a timer for each kernel
	usetimer timer_heatbath_odd;
	usetimer timer_heatbath_even;
	usetimer timer_overrelax_odd;
	usetimer timer_overrelax_even;

	/**
	 * Return the timer connected to a specific kernel.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual usetimer* get_timer(const char * in) override;

	/**
	 * Print the profiling information to a file.
	 *
	 * @param filename Name of file where data is appended.
	 * @param parameters inputparameters
	 */
	void virtual print_profiling(std::string filename, int number) override;
#endif

	/**
	 * Return amount of bytes read and written by a specific kernel per call.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual size_t get_read_write_size(const char * in) override;

	/**
	 * Return amount of Floating point operations performed by a specific kernel per call.
	 * NOTE: this is meant to be the "netto" amount in order to be comparable.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual uint64_t get_flop_size(const char * in) override;


};

#endif //OPENCLMODULEHEATBATHH
