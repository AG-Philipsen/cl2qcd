/** @file
 * Basic OpenCL functionality
 */
#ifndef _OPENCLMODULEMOLECULARDYNAMICSH_
#define _OPENCLMODULEMOLECULARDYNAMICSH_

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
#include "types_hmc.h"

#include "exceptions.h"

#include "meta/counter.hpp"
#include "hardware/buffers/plain.hpp"
#include "hardware/buffers/prng_buffer.hpp"
#include "hardware/buffers/su3.hpp"
#include "hardware/buffers/spinor.hpp"
#include "hardware/buffers/gaugemomentum.hpp"

/**
 * An OpenCL device
 *
 * This class wraps all operations on a device. Operations are always specific, e.g. each kernel and copy operation
 * has it's own wrapper function.
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Opencl_Module_Molecular_Dynamics : public Opencl_Module {
public:
	friend hardware::Device;

	virtual ~Opencl_Module_Molecular_Dynamics();

	///////////////////////////////////////////////////
	//Methods on device
	void md_update_gaugemomentum_device(const hardware::buffers::Gaugemomentum *, const hardware::buffers::Gaugemomentum *, hmc_float eps);
	void md_update_gaugefield_device(hmc_float eps);
	void md_update_gaugefield_device(const hardware::buffers::Gaugemomentum *, const hardware::buffers::SU3 *, hmc_float eps);
	void gauge_force_device(const hardware::buffers::SU3 * gf, const hardware::buffers::Gaugemomentum * out);
	void gauge_force_tlsym_device(const hardware::buffers::SU3 * gf, const hardware::buffers::Gaugemomentum * out);
	void fermion_force_device(const hardware::buffers::Plain<spinor> * Y, const hardware::buffers::Plain<spinor> * X, hmc_float kappa = ARG_DEF);
	void fermion_force_device(const hardware::buffers::Plain<spinor> * Y, const hardware::buffers::Plain<spinor> * X, const hardware::buffers::SU3 *, const hardware::buffers::Gaugemomentum *, hmc_float kappa = ARG_DEF);
	void fermion_force_eo_device(const hardware::buffers::Spinor * Y, const hardware::buffers::Spinor * X, int evenodd, hmc_float kappa = ARG_DEF);
	void fermion_force_eo_device(const hardware::buffers::Spinor * Y, const hardware::buffers::Spinor * X, const hardware::buffers::SU3 *, const hardware::buffers::Gaugemomentum *, int evenodd, hmc_float kappa = ARG_DEF);
	void stout_smeared_fermion_force_device(std::vector<const hardware::buffers::SU3 *>& gf_intermediate);

protected:

	/**
	 * comutes work-sizes for a kernel
	 * @todo autotune
	 * @param ls local-work-size
	 * @param gs global-work-size
	 * @param num_groups number of work groups
	 * @param name name of the kernel for possible autotune-usage, not yet used!!
	 */
	virtual void get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const override;

	/**
	 * Print the profiling information to a file.
	 *
	 * @param filename Name of file where data is appended.
	 */
	void virtual print_profiling(const std::string& filename, int number) const override;

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
	/**
	 * Constructor.
	 *
	 * @param[in] params points to an instance of inputparameters
	 */
	Opencl_Module_Molecular_Dynamics(const meta::Inputparameters& params, hardware::Device * device);

	/**
	 * Collect the kernels for OpenCL.
	 */
	void fill_kernels();
	/**
	 * Clear out the kernels,
	 */
	void clear_kernels();

	ClSourcePackage basic_hmc_code;

	//kernels
	cl_kernel md_update_gaugefield;
	cl_kernel md_update_gaugemomenta;
	cl_kernel gauge_force;
	cl_kernel gauge_force_tlsym;
	cl_kernel fermion_force;
	cl_kernel fermion_force_eo;
	cl_kernel stout_smear_fermion_force;
};

#endif //OPENCLMODULEMOLECULARDYNAMICSH
