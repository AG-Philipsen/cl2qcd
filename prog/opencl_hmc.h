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
	 * Collect a vector of kernel file names.
	 * Virtual method, allows to include more kernel files in inherited classes.
	 */
	virtual hmc_error fill_kernels_file ();
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
	virtual hmc_error fill_kernels();


	/**
	 * Initialize the OpenCL device including fermion capabilities
	 *
	 * @param wanted The OpenCL device type to be used, e.g. CL_DEVICE_TYPE_CPU or CL_DEVICE_TYPE_GPU
	 * @param timer The timer to use for reporting execution time
	 * @param parameters The parsed input parameters
	 * @return Error code as defined in hmcerrs.h:
	 *         @li HMC_OCLERROR if OpenCL initialization / operations fail
	 *         @li HMC_FILEERROR if one of the kernel files cannot be opened
	 *         @li HMC_SUCCESS otherwise
	 */
	virtual hmc_error init(cl_device_type wanted_device_type, usetimer* timer, inputparameters* parameters);

	hmc_error finalize_hmc();

	hmc_error init_hmc_variables(inputparameters* parameters, usetimer* timer);
	
	////////////////////////////////////////////////////
	//Methods needed for the HMC-algorithm
	
	hmc_error generate_gaussian_gaugemomenta_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer);
	hmc_error generate_gaussian_spinorfield_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer);
	hmc_error md_update_spinorfield_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer);;
	hmc_error leapfrog_device(hmc_float tau, int steps1, int steps2, const size_t local_work_size, const size_t global_work_size, usetimer * timer);
	hmc_error force_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer);
	hmc_error hamiltonian_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer);
	hmc_error calc_spinorfield_init_energy_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer);
	hmc_error md_update_gaugemomentum_device(hmc_float eps, const size_t local_work_size, const size_t global_work_size, usetimer * timer);
	hmc_error md_update_gaugefield_device(hmc_float eps, const size_t local_work_size, const size_t global_work_size, usetimer * timer);
	hmc_error set_zero_clmem_force_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer);
	hmc_error gauge_force_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer);
	hmc_error fermion_force_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer);
	
	
	
	////////////////////////////////////////////////////
	//copying
	//Methods to copy new and old fields... these can be optimized!!
	hmc_error copy_gaugefield_old_new_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer);
	hmc_error copy_gaugemomenta_old_new_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer);
	hmc_error copy_gaugefield_new_old_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer);
	hmc_error copy_gaugemomenta_new_old_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer);
	
	hmc_error get_deltah_from_device(hmc_float * out, const size_t local_work_size, const size_t global_work_size, usetimer * timer);
	private:
		//kernels
		cl_kernel generate_gaussian_spinorfield;
		cl_kernel generate_gaussian_gaugemomenta;
		cl_kernel md_update_gaugefield;
		cl_kernel md_update_gaugemomenta;
		cl_kernel gauge_force;
		cl_kernel fermion_force;
		cl_kernel s_gauge;
		cl_kernel s_fermion;
		cl_kernel set_zero_gaugemomentum;
		
		//variables
		//initial energy of the (gaussian) spinorfield
		cl_mem clmem_energy_init;
		//DeltaH
		cl_mem clmem_deltah;
		//new and old gaugemomentum, new gaugefield
		cl_mem clmem_p;
		cl_mem clmem_new_p;
		cl_mem clmem_new_u;
		//force field
		cl_mem clmem_force;
		//inverted spinorfield
		cl_mem clmem_phi_inv;
		
};
#endif // _MYOPENCLHMCH_
