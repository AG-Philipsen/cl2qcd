/** @file
 * OpenCL device including heatbath.
 */
#ifndef _MYOPENCLHEATBATHH_
#define _MYOPENCLHEATBATHH_

#include "opencl.h"


/**
 * An OpenCL device for heatbath.
 *
 * This class wraps all operations on a device. Operations are always specific, e.g. each kernel and copy operation
 * has it's own wrapper function.
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Opencl_heatbath : public Opencl {
public:

	//Calculations, calls to kernels
	/**
	 * Perform one heatbath step.
	 */
	hmc_error run_heatbath(const hmc_float beta);

	/**
	 * Perform one overrelaxation step.
	 */
	hmc_error run_overrelax(const hmc_float beta);

	/**
	 * Calculate plaquette and polyakov.
	 *
	 * @param[out] plaq Storage for result of plaquette calculation
	 * @param[out] tplaq Storage for result of plaquette calculation
	 * @param[out] splaq Storage for result of plaquette calculation
	 * @param[out] pol Storage for result of polyakov calculation
	 * @param[in,out] timer1 Timer into which to aggregate plaquette calculation time
	 * @param[in,out] timer2 Timer into which to aggregate polyakov calculation time
	 * @return Error code as defined in hmcerrs.h:
	 *         @li HMC_OCLERROR if OpenCL operations fail
	 *         @li HMC_SUCCESS otherwise
	 */
	//protected:

	/**
	 * Collect the compiler options for OpenCL.
	 * Virtual method, allows to include more options in inherited classes.
	 */
	virtual hmc_error fill_collect_options(stringstream* collect_options);

	/**
	 * Collect the buffers to generate for OpenCL.
	 * Virtual method, allows to include more buffers in inherited classes.
	 */
	virtual hmc_error fill_buffers();

	/**
	 * Collect the kernels for OpenCL.
	 * Virtual method, allows to include more kernels in inherited classes.
	 */
	virtual void fill_kernels();

	//private:


	/**
	 * Clear out the kernels,
	 * Virtual method, allows to clear additional kernels in inherited classes.
	 */
	virtual hmc_error clear_kernels();

	/**
	 * Clear out the buffers,
	 * Virtual method, allows to clear additional buffers in inherited classes.
	 */
	virtual hmc_error clear_buffers();

	///////////////////////////////////////////////////////////
	//LZ what follows should eventually be private
	//heatbath variables



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


protected:

private:

};

#endif /* _MYOPENCLHEATBATHH_ */
