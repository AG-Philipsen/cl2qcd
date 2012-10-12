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
	Opencl_Module_Hmc(const meta::Inputparameters& params, hardware::Device * device, meta::Counter * inversions0, meta::Counter * inversions1,
	                  meta::Counter * inversions_mp0, meta::Counter * inversions_mp1)
		: Opencl_Module_Fermions(params, device), new_u(get_gaugefield()->get_elements(), device), gaugemomentum_buf_size(0),
		  inversions0(inversions0), inversions1(inversions1), inversions_mp0(inversions_mp0), inversions_mp1(inversions_mp1),
		  clmem_phi_inv(meta::get_spinorfieldsize(params), device),
		  clmem_phi_inv_eo(meta::get_eoprec_spinorfieldsize(params), device),
		  clmem_phi(meta::get_spinorfieldsize(params), device),
		  clmem_phi_mp(meta::get_spinorfieldsize(params), device),
		  clmem_phi_eo(meta::get_eoprec_spinorfieldsize(params), device),
		  clmem_phi_mp_eo(meta::get_eoprec_spinorfieldsize(params), device)
		{ };

	// OpenCL specific methods needed for building/compiling the OpenCL program
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
	 * comutes work-sizes for a kernel
	 * @todo autotune
	 * @param ls local-work-size
	 * @param gs global-work-size
	 * @param num_groups number of work groups
	 * @param name name of the kernel for possible autotune-usage, not yet used!!
	 */
	virtual void get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const override;

	////////////////////////////////////////////////////
	//get members
	cl_mem get_clmem_p();
	cl_mem get_clmem_new_p();
	const hardware::buffers::SU3 * get_new_u();
	const hardware::buffers::ScalarBuffer<spinor> * get_clmem_phi();
	const hardware::buffers::Spinor * get_clmem_phi_eo();
	const hardware::buffers::ScalarBuffer<spinor> * get_clmem_phi_mp();
	const hardware::buffers::Spinor * get_clmem_phi_mp_eo();
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
	void md_update_gaugefield_device(cl_mem, const hardware::buffers::SU3 *, hmc_float eps);
	void set_zero_clmem_force_device();
	void gauge_force_device();
	void gauge_force_device(const hardware::buffers::SU3 * gf, cl_mem out);
	void gauge_force_tlsym_device();
	void gauge_force_tlsym_device(const hardware::buffers::SU3 * gf, cl_mem out);
	void fermion_force_device(const hardware::buffers::ScalarBuffer<spinor> * Y, const hardware::buffers::ScalarBuffer<spinor> * X, hmc_float kappa = ARG_DEF);
	void fermion_force_device(const hardware::buffers::ScalarBuffer<spinor> * Y, const hardware::buffers::ScalarBuffer<spinor> * X, const hardware::buffers::SU3 *, cl_mem, hmc_float kappa = ARG_DEF);
	void fermion_force_eo_device(const hardware::buffers::Spinor * Y, const hardware::buffers::Spinor * X, int evenodd, hmc_float kappa = ARG_DEF);
	void fermion_force_eo_device(const hardware::buffers::Spinor * Y, const hardware::buffers::Spinor * X, const hardware::buffers::SU3 *, cl_mem, int evenodd, hmc_float kappa = ARG_DEF);
	void stout_smeared_fermion_force_device(std::vector<const hardware::buffers::SU3 *>& gf_intermediate);
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

	/**
	 * Print the profiling information to a file.
	 *
	 * @param filename Name of file where data is appended.
	 */
	void virtual print_profiling(const std::string& filename, int number) override;

	/**
	 * Return amount of bytes read and written by a specific kernel per call.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual size_t get_read_write_size(const std::string& in) const override;

	/**
	 * Return amount of Floating point operations performed by a specific kernel per call.
	 * NOTE: this is meant to be the "netto" amount in order to be comparable.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual uint64_t get_flop_size(const std::string& in) const override;

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
	const hardware::buffers::SU3 new_u;
	//force field
	cl_mem clmem_force;
	//inverted spinorfield
	const hardware::buffers::ScalarBuffer<spinor> clmem_phi_inv;
	const hardware::buffers::Spinor clmem_phi_inv_eo;
	//D(gaussian spinorfield)
	const hardware::buffers::ScalarBuffer<spinor> clmem_phi;
	const hardware::buffers::ScalarBuffer<spinor> clmem_phi_mp;
	const hardware::buffers::Spinor clmem_phi_eo;
	const hardware::buffers::Spinor clmem_phi_mp_eo;

	ClSourcePackage basic_hmc_code;

	size_t gaugemomentum_buf_size;

	meta::Counter * const inversions0;
	meta::Counter * const inversions1;
	meta::Counter * const inversions_mp0;
	meta::Counter * const inversions_mp1;
};

#endif //OPENCLMODULEHMCH
