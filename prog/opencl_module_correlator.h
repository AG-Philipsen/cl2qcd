/** @file
 * Basic OpenCL functionality
 */
#ifndef _OPENCLMODULECORRELATORH_
#define _OPENCLMODULECORRELATORH_

#include "opencl_module.h"

#include "hardware/buffers/plain.hpp"
#include "hardware/buffers/prng_buffer.hpp"
#include "types_fermions.h"

/**
 * An OpenCL device
 *
 * This class wraps all operations on a device. Operations are always specific, e.g. each kernel and copy operation
 * has it's own wrapper function.
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Opencl_Module_Correlator : public Opencl_Module {
public:
	friend hardware::Device;

	virtual ~Opencl_Module_Correlator();

	void create_point_source_device(const hardware::buffers::Plain<spinor> * inout, int i, int spacepos, int timepos);

	void create_volume_source_device(const hardware::buffers::Plain<spinor> * inout, const hardware::buffers::PRNGBuffer * prng);

	void create_timeslice_source_device(const hardware::buffers::Plain<spinor> * inout, const hardware::buffers::PRNGBuffer * prng, const int timeslice);

	void create_zslice_source_device(const hardware::buffers::Plain<spinor> * inout, const hardware::buffers::PRNGBuffer * prng, const int zslice);

	/**
	 * Calculate specific correlator on device.
	 * This function is overloaded depending on whether one needs the source for the calculation or not.
	 */
	void correlator_device(const cl_kernel correlator_kernel, const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<hmc_float> * correlator);
	void correlator_device(const cl_kernel correlator_kernel, const hardware::buffers::Plain<spinor> * in, const hardware::buffers::Plain<spinor> * source, const hardware::buffers::Plain<hmc_float> * correlator);

	/**
	 * Get kernel for correlator indicated by which
	 * @param[in] which string that identifies the correlator (ps or sc, vx, vy, vz, ax, ay, az)
	 * @return correlator_kernel
	 */
	cl_kernel get_correlator_kernel(std::string which);

	/////////////////////////////////////////////////
	//functions to get private variables
	cl_mem get_clmem_corr();

	/**
	 * Print the profiling information to a file.
	 *
	 * @param filename Name of file where data is appended.
	 * @param number task-id
	 */
	void virtual print_profiling(const std::string& filename, int number) const override;

protected:
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

	/**
	 * comutes work-sizes for a kernel
	 * @todo autotune
	 * @param ls local-work-size
	 * @param gs global-work-size
	 * @param num_groups number of work groups
	 * @param name name of the kernel for possible autotune-usage, not yet used!!
	 */
	virtual void get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const override;

private:
	/**
	 * Default constructor, does nothing but make sure some pointer point to 0.
	 *
	 */
	Opencl_Module_Correlator(const meta::Inputparameters& params, hardware::Device * device);

	/**
	 * Collect the kernels for OpenCL.
	 */
	void fill_kernels();
	/**
	 * Clear out the kernels,
	 */
	void clear_kernels();

	////////////////////////////////////
	//kernels

	cl_kernel create_point_source;
	cl_kernel create_volume_source;
	cl_kernel create_timeslice_source;
	cl_kernel create_zslice_source;

	//Observables
	//scalar correlators
	cl_kernel correlator_ps;
	cl_kernel correlator_sc;
	//vector correlators, physical basis: Gamma_4 * Gamma_{x,y,z}
	cl_kernel correlator_vx;
	cl_kernel correlator_vy;
	cl_kernel correlator_vz;
	//axial vector correlators, physical basis: Gamma_5 * Gamma_4 * Gamma_{x,y,z}
	cl_kernel correlator_ax;
	cl_kernel correlator_ay;
	cl_kernel correlator_az;

	ClSourcePackage basic_correlator_code;
};

#endif //OPENCLMODULECORRELATORH
