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
#include "host_operations_gaugefield.h"
#include "globaldefs.h"
#include "types.h"
#include "types_fermions.h"
#include "host_use_timer.h"
#include "opencl_compiler.hpp"

#include "opencl_module.h"
#include "opencl_module_ran.h"
#include "opencl_module_spinors.h"
#include "opencl_module_fermions.h"
#include "types_hmc.h"

#include "exceptions.h"

#include "meta/counter.hpp"

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

	/**
	 * Empty constructor.
	 *
	 * @param[in] params points to an instance of inputparameters
	 */
	Opencl_Module_Hmc(const meta::Inputparameters& params, meta::Counter * inversions0, meta::Counter * inversions1,
	                  meta::Counter * inversions_mp0, meta::Counter * inversions_mp1)
		: Opencl_Module_Fermions(params), gaugemomentum_buf_size(0),
		  inversions0(inversions0), inversions1(inversions1), inversions_mp0(inversions_mp0), inversions_mp1(inversions_mp1) { }

	// OpenCL specific methods needed for building/compiling the OpenCL program
	/**
	 * Collect the compiler options for OpenCL.
	 * Virtual method, allows to include more options in inherited classes.
	 */
	virtual void fill_collect_options(std::stringstream* collect_options);
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
	cl_mem get_clmem_phi_eo();
	cl_mem get_clmem_phi_mp();
	cl_mem get_clmem_phi_mp_eo();
	cl_mem get_clmem_s_fermion_init();
  	cl_mem get_clmem_s_fermion_mp_init();

	////////////////////////////////////////////////////
	//Methods needed for the HMC-algorithm
	void md_update_spinorfield(hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void md_update_spinorfield_mp(usetimer * solvertimer);
	void generate_spinorfield_gaussian();
	hmc_observables metropolis(hmc_float rnd, hmc_float beta);
	void calc_spinorfield_init_energy(cl_mem dest);
	void calc_gauge_force();
	void calc_fermion_force(usetimer * solvertimer, hmc_float kappa = ARG_DEF, hmc_float mubar = ARG_DEF);
	void calc_fermion_force_detratio(usetimer * solvertimer);

	///////////////////////////////////////////////////
	//Methods on device
	void set_float_to_gaugemomentum_squarenorm_device(cl_mem in, cl_mem out);
	void generate_gaussian_gaugemomenta_device();
	void generate_gaussian_spinorfield_device();
	void generate_gaussian_spinorfield_eo_device();
	void md_update_gaugemomentum_device(hmc_float eps);
	void md_update_gaugemomentum_device(cl_mem, cl_mem, hmc_float eps);
	void md_update_gaugefield_device(hmc_float eps);
	void md_update_gaugefield_device(cl_mem, cl_mem, hmc_float eps);
	void set_zero_clmem_force_device();
	void gauge_force_device();
	void gauge_force_device(cl_mem gf, cl_mem out);
	void gauge_force_tlsym_device();
	void gauge_force_tlsym_device(cl_mem gf, cl_mem out);
	void fermion_force_device(cl_mem Y, cl_mem X, hmc_float kappa = ARG_DEF);
	void fermion_force_eo_device(cl_mem Y, cl_mem X, int evenodd, hmc_float kappa = ARG_DEF);
	void stout_smeared_fermion_force_device(cl_mem * gf_intermediate);
	hmc_float calc_s_fermion();
	hmc_float calc_s_fermion_mp();

	/**
	 * Query the size required for storage of a gaugemomentum buffer
	 * in the format used by the implementation.
	 *
	 * @return The buffer size in bytes
	 */
	size_t get_gaugemomentum_buffer_size();

	/**
	 * Import data from the gaugemomenta array into the given buffer.
	 *
	 * The data in the buffer will be stored in the device specific format.
	 *
	 * @param[out] dest The buffer to write to in the device specific format
	 * @param[in] data The data to write to the buffer
	 */
	void importGaugemomentumBuffer(const cl_mem dest, const ae * const data);
	/**
	 * Export data from the given buffer into a normal gaugemomentum array.
	 *
	 * The data in the buffer is assumed to be in the device specific format.
	 *
	 * @param[out] dest An array that the buffer data can be written to.
	 * @param[in] data A buffer containing the data in the device specific format.
	 */
	void exportGaugemomentumBuffer(ae * const dest, const cl_mem buf);

protected:

#ifdef _PROFILING_

	usetimer timer_generate_gaussian_spinorfield;
	usetimer timer_generate_gaussian_spinorfield_eo;
	usetimer timer_generate_gaussian_gaugemomenta;
	usetimer timer_md_update_gaugefield;
	usetimer timer_md_update_gaugemomenta;
	usetimer timer_gauge_force;
	usetimer timer_gauge_force_tlsym;
	usetimer timer_fermion_force;
	usetimer timer_fermion_force_eo;
	usetimer timer_set_zero_gaugemomentum;
	usetimer timer_gaugemomentum_squarenorm;
	usetimer timer_stout_smear_fermion_force;

	/**
	 * Return the timer connected to a specific kernel.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual usetimer* get_timer(const char * in);

	/**
	 * Print the profiling information to a file.
	 *
	 * @param filename Name of file where data is appended.
	 */
	void virtual print_profiling(std::string filename, int number);

#endif

	/**
	 * Return amount of bytes read and written by a specific kernel per call.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual size_t get_read_write_size(char * in);

	/**
	 * Return amount of Floating point operations performed by a specific kernel per call.
	 * NOTE: this is meant to be the "netto" amount in order to be comparable.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual uint64_t get_flop_size(const char * in);

private:

	//kernels
	cl_kernel generate_gaussian_spinorfield;
	cl_kernel generate_gaussian_spinorfield_eo;
	cl_kernel generate_gaussian_gaugemomenta;
	cl_kernel md_update_gaugefield;
	cl_kernel md_update_gaugemomenta;
	cl_kernel gauge_force;
	cl_kernel gauge_force_tlsym;
	cl_kernel fermion_force;
	cl_kernel fermion_force_eo;
	cl_kernel stout_smear_fermion_force;
	cl_kernel set_zero_gaugemomentum;
	cl_kernel gaugemomentum_squarenorm;
	cl_kernel gaugemomentum_convert_to_soa;
	cl_kernel gaugemomentum_convert_from_soa;

	//variables
	//initial energy of the (gaussian) spinorfield
	cl_mem clmem_s_fermion_init;
	cl_mem clmem_s_fermion_mp_init;
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
	cl_mem clmem_phi_inv_eo;
	//D(gaussian spinorfield)
	cl_mem clmem_phi;
	cl_mem clmem_phi_mp;
	cl_mem clmem_phi_eo;
	cl_mem clmem_phi_mp_eo;

	ClSourcePackage basic_hmc_code;

	size_t gaugemomentum_buf_size;

	meta::Counter * const inversions0;
	meta::Counter * const inversions1;
	meta::Counter * const inversions_mp0;
	meta::Counter * const inversions_mp1;
};

#endif //OPENCLMODULEHMCH
