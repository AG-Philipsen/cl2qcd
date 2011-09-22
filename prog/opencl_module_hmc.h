/** @file
 * Basic OpenCL functionality
 */
#ifndef _OPENCLMODULEHMCH_
#define _OPENCLMODULEHMCH_

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
#include "types_fermions.h"
#include "host_use_timer.h"
#include "inputparameters.h"
#include "opencl_compiler.hpp"

#include "opencl_module.h"
#include "opencl_module_ran.h"
#include "opencl_module_spinors.h"
#include "opencl_module_fermions.h"
#include "types_hmc.h"

#include "exceptions.h"

/**
 * An OpenCL device
 *
 * This class wraps all operations on a device. Operations are always specific, e.g. each kernel and copy operation
 * has it's own wrapper function.
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Opencl_Module_Hmc : public Opencl_Module_Fermions {
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

	////////////////////////////////////////////////////
	//get members
	cl_mem get_clmem_p();
	cl_mem get_clmem_new_p();
	cl_mem get_clmem_new_u();
	cl_mem get_clmem_phi();
	cl_mem get_clmem_phi_eoprec();
	
	////////////////////////////////////////////////////
	//Methods needed for the HMC-algorithm
	void md_update_spinorfield();
	hmc_observables metropolis(hmc_float rnd, hmc_float beta);
	void calc_spinorfield_init_energy();
	void calc_gauge_force();
	void calc_fermion_force(usetimer * solvertimer);

	///////////////////////////////////////////////////
	//Methods on device
	void set_float_to_gaugemomentum_squarenorm_device(cl_mem in, cl_mem out);
	void generate_gaussian_gaugemomenta_device();
	void generate_gaussian_spinorfield_device();
	void generate_gaussian_spinorfield_eoprec_device();
	void md_update_gaugemomentum_device(hmc_float eps);
	void md_update_gaugefield_device(hmc_float eps);
	void set_zero_clmem_force_device();
	void gauge_force_device();
	void fermion_force_device();
	void fermion_force_eoprec_device(cl_mem Y, cl_mem X);
	void stout_smeared_fermion_force_device();
	
protected:

#ifdef _PROFILING_

	usetimer timer_generate_gaussian_spinorfield;
	usetimer timer_generate_gaussian_spinorfield_eoprec;
	usetimer timer_generate_gaussian_gaugemomenta;
	usetimer timer_md_update_gaugefield;
	usetimer timer_md_update_gaugemomenta;
	usetimer timer_gauge_force;
	usetimer timer_fermion_force;
	usetimer timer_fermion_force_eoprec;
	usetimer timer_set_zero_gaugemomentum;
	usetimer timer_gaugemomentum_squarenorm;
	usetimer timer_stout_smear_fermion_force;
	
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
	 */
	void virtual print_profiling(std::string filename);

#endif

protected:
private:

	//kernels
	cl_kernel generate_gaussian_spinorfield;
	cl_kernel generate_gaussian_spinorfield_eoprec;
	cl_kernel generate_gaussian_gaugemomenta;
	cl_kernel md_update_gaugefield;
	cl_kernel md_update_gaugemomenta;
	cl_kernel gauge_force;
	cl_kernel fermion_force;
	cl_kernel fermion_force_eoprec;
	cl_kernel stout_smear_fermion_force;
	cl_kernel set_zero_gaugemomentum;
	cl_kernel gaugemomentum_squarenorm;

	//variables
	//initial energy of the (gaussian) spinorfield
	cl_mem clmem_energy_init;
	//squarenorm temps
	cl_mem clmem_p2;
	cl_mem clmem_new_p2;
	cl_mem clmem_s_fermion;
	//new and old gaugemomentum, new gaugefield
	cl_mem clmem_p;
	cl_mem clmem_new_p;
	cl_mem clmem_new_u;
	//force field
	cl_mem clmem_force;
	//inverted spinorfield
	cl_mem clmem_phi_inv;
	cl_mem clmem_phi_inv_eoprec;
	//D(gaussian spinorfield)
	cl_mem clmem_phi;
	cl_mem clmem_phi_eoprec;
	//these are the odd vectors for the force-calculation
	cl_mem clmem_x_odd;
	cl_mem clmem_y_odd;

	ClSourcePackage basic_hmc_code;

};

#endif //OPENCLMODULEHMCH
