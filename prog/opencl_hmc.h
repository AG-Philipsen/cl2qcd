/** @file
 * OpenCL device managment and everything executed on them -- including HMC-algorithm.
 */
#ifndef _MYOPENCLHMCH_
#define _MYOPENCLHMCH_

//these are the headers of the mother-classes
#include "opencl.h"
#include "opencl_fermions.h"
//CP: this includes the struct-definitions for the spinors...
#include "types_hmc.h"

/**
 * An OpenCL device for the HMC algorithm.
 *
 * This class wraps all operations on a device. Inherited from classes Opencl and Opencl_fermions.
 */
class Opencl_hmc : public Opencl_fermions {
public:
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
	/**
	 * Initialize the OpenCL device including fermion capabilities
	 *
	 * @param wanted The OpenCL device type to be used, e.g. CL_DEVICE_TYPE_CPU or CL_DEVICE_TYPE_GPU
	 * @param parameters The parsed input parameters
	 * @param nstates Number of random states
	 * @return Error code as defined in hmcerrs.h:
	 *         @li HMC_OCLERROR if OpenCL initialization / operations fail
	 *         @li HMC_FILEERROR if one of the kernel files cannot be opened
	 *         @li HMC_SUCCESS otherwise
	 */
	virtual hmc_error init(cl_device_type wanted_device_type, inputparameters* parameters, int nstates);

	hmc_error finalize_hmc();

	////////////////////////////////////////////////////
	//Methods needed for the HMC-algorithm
	void generate_gaussian_gaugemomenta_device();
	void generate_gaussian_spinorfield_device();
	void md_update_spinorfield_device();
	void leapfrog_device(hmc_float tau, int steps1, int steps2, usetimer *copy_to, usetimer * copy_on, usetimer * solvertimer);
	hmc_observables metropolis(hmc_float rnd, hmc_float beta, usetimer * timer);
	void calc_spinorfield_init_energy_device();
	void md_update_gaugemomentum_device(hmc_float eps);
	void md_update_gaugefield_device(hmc_float eps);
	void set_zero_clmem_force_device();
	void gauge_force_device();
	void fermion_force_device();
	void stout_smeared_fermion_force_device();
	void set_float_to_gaugemomentum_squarenorm_device(cl_mem in, cl_mem out);
	void calc_total_force(usetimer *copy_to, usetimer * copy_on, usetimer * solvertimer);
	void calc_gauge_force();
	void calc_fermion_force(usetimer *copy_to, usetimer * copy_on, usetimer * solvertimer);
	
	////////////////////////////////////////////////////
	//get members
	cl_mem get_clmem_p();
	cl_mem get_clmem_new_p();
	cl_mem get_clmem_new_u();
	cl_mem get_clmem_phi();
	
#ifdef _PROFILING_
	//CP: if PROFILING is activated, one needs a timer for each kernel
	//fermionmatrix
	usetimer timer_generate_gaussian_spinorfield;
	usetimer timer_generate_gaussian_gaugemomenta;
	usetimer timer_md_update_gaugefield;
	usetimer timer_md_update_gaugemomenta;
	usetimer timer_gauge_force;
	usetimer timer_fermion_force;
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
	 * @param parameters inputparameters
	 */
	void virtual print_profiling(std::string filename);	
#endif		
	
private:
	//kernels
	cl_kernel generate_gaussian_spinorfield;
	cl_kernel generate_gaussian_gaugemomenta;
	cl_kernel md_update_gaugefield;
	cl_kernel md_update_gaugemomenta;
	cl_kernel gauge_force;
	cl_kernel fermion_force;
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
	//D(gaussian spinorfield)
	cl_mem clmem_phi;

	ClSourcePackage basic_hmc_code;

};
#endif // _MYOPENCLHMCH_
