/** @file
 * Opencl code working on the gaugefield.
 */
#ifndef _OPENCL_MODULE_GAUGEFIELD_H_
#define _OPENCL_MODULE_GAUGEFIELD_H_

#include "opencl_module.h"
#include "types.h"
#include "meta/inputparameters.hpp"
#include "hardware/device.hpp"
#include "opencl_compiler.hpp"
#include "hardware/buffers/su3.hpp"
#include "hardware/buffers/plain.hpp"

/**
 * An OpenCL device
 *
 * This class wraps all gaugefield operations on a device. Each kernel
 * has it's own wrapper function.
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Opencl_Module_Gaugefield : public Opencl_Module {

public:
	/**
	 * Empty constructor.
	 *
	 * @param[in] params points to an instance of inputparameters
	 */
	Opencl_Module_Gaugefield(const meta::Inputparameters& params, hardware::Device * device);
	/**
	 * Destructor.
	 */
	virtual ~Opencl_Module_Gaugefield();

	/**
	 * Get a pointer to the gaugefield buffer
	 * @return ocl_gaugefield OpenCL buffer with gaugefield
	 */
	const hardware::buffers::SU3 * get_gaugefield();

	// methods which actually calculate something
	/**
	 * Calculate plaquette and polyakov of a specific gaugefield (on device).
	 *
	 * @param[out] plaq Storage for result of plaquette calculation
	 * @param[out] tplaq Storage for result of plaquette calculation
	 * @param[out] splaq Storage for result of plaquette calculation
	 * @param[out] pol Storage for result of polyakov calculation
	 */
	void gaugeobservables(hmc_float * const plaq, hmc_float * const tplaq, hmc_float * const splaq, hmc_complex * const pol);
	/**
	 * Calculate plaquette and polyakov of a specific gaugefield (on device).
	 *
	 * @param[in]  gf    The gaugefield on which to compute the observables
	 * @param[out] plaq  Storage for result of plaquette calculation
	 * @param[out] tplaq Storage for result of plaquette calculation
	 * @param[out] splaq Storage for result of plaquette calculation
	 * @param[out] pol   Storage for result of polyakov calculation
	 *
	 * @todo Should not be public
	 */
	void gaugeobservables(const hardware::buffers::SU3 * gf, hmc_float * const plaq, hmc_float * const tplaq, hmc_float * const splaq, hmc_complex * const pol);
	/**
	 * Calculate rectangles of a specific gaugefield (on device).
	 *
	 * @param[in] gf gaugefield to measure on
	 * @param[out] plaq Storage for result of rectangles calculation
	 */
	void gaugeobservables_rectangles(const hardware::buffers::SU3 * gf, hmc_float * const rect);
	/**
	 * Calculate plaquette for a specific gaugefield (on device).
	 *
	 * @param[in] gf gaugefield to measure on
	 */
	void plaquette_device(const hardware::buffers::SU3 * gf, const hardware::buffers::Plain<hmc_float> * plaq, const hardware::buffers::Plain<hmc_float> * tplaq, const hardware::buffers::Plain<hmc_float> * splaq);
	/**
	 * Calculate rectangles for a specific gaugefield (on device).
	 *
	 * @param[in] gf gaugefield to measure on
	 */
	void rectangles_device(const hardware::buffers::SU3 * gf, const hardware::buffers::Plain<hmc_float> *);
	/**
	 * Calculate Polyakov loop for a specific gaugefield (on device).
	 *
	 * @param[in] gf gaugefield to measure on
	 */
	void polyakov_device(const hardware::buffers::SU3 * gf, const hardware::buffers::Plain<hmc_complex> *);

	/**
	 * Print the profiling information to a file.
	 *
	 * @param filename Name of file where data is appended.
	 * @param number task-id
	 */
	void virtual print_profiling(const std::string& filename, int number) const override;

	/**
	 * This applies stout smearing to a gaugefield
	 */
	void smear_gaugefield(const hardware::buffers::SU3 * gf, const std::vector<const hardware::buffers::SU3 *>& gf_intermediate);
	void stout_smear_device(const hardware::buffers::SU3 * in, const hardware::buffers::SU3 * out);

	/**
	 * This replaces the stout smeared gaugefield with the unsmeared one
	 */
	void unsmear_gaugefield(const hardware::buffers::SU3 * gf);

	/**
	 * Import the gaugefield data into the OpenCL buffer using the device
	 * specific storage format.
	 *
	 * @param[in]  data       The gaugefield data to import into the OpenCL buffer
	 * @param[out] gaugefield The OpenCL buffer to writ the gaugefield data to in the device specific format
	 */
	void importGaugefield(const Matrixsu3 * const data);

	/**
	 * Import the gaugefield data into the OpenCL buffer using the device
	 * specific storage format.
	 *
	 * @param[out] gaugefield The OpenCL buffer to writ the gaugefield data to in the device specific format
	 * @param[in]  data       The gaugefield data to import into the OpenCL buffer
	 *
	 * @todo should not be public
	 */
	void importGaugefield(const hardware::buffers::SU3 * gaugefield, const Matrixsu3 * const data);

	/**
	 * Export the gaugefield from the OpenCL buffer, that uses a device
	 * specific storage format, into the given pointer using the generic
	 * storage format.
	 *
	 * @param[out] dest The array to store the gaugefield in
	 */
	void exportGaugefield(Matrixsu3 * const dest);

protected:
	/**
	 * comutes work-sizes for a kernel
	 * @todo autotune
	 * @param ls local-work-size
	 * @param gs global-work-size
	 * @param num_groups number of work groups
	 * @param dev_type type of device on which the kernel should be executed
	 * @param name name of the kernel for possible autotune-usage, not yet used!!
	 */
	virtual void get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const override;

	/**
	 * Return amount of Floating point operations performed by a specific kernel per call.
	 * NOTE: this is meant to be the "netto" amount in order to be comparable.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual uint64_t get_flop_size(const std::string& in) const override;

	/**
	 * Return amount of bytes read and written by a specific kernel per call.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual size_t get_read_write_size(const std::string& in) const override;

	/**
	 * A set of source files used by all kernels.
	 */
	ClSourcePackage basic_opencl_code;

private:

	const hardware::buffers::SU3 gaugefield;

	//this is used to save the unsmeared gaugefield if smearing is used
	const hardware::buffers::SU3 gf_unsmeared;

	//since this is only applicated to the gaugefield, this should be here...
	cl_kernel stout_smear;

	cl_kernel plaquette;
	cl_kernel plaquette_reduction;
	cl_kernel rectangles;
	cl_kernel rectangles_reduction;
	cl_kernel polyakov;
	cl_kernel polyakov_reduction;
	cl_kernel convertGaugefieldToSOA;
	cl_kernel convertGaugefieldFromSOA;

	void convertGaugefieldToSOA_device(const hardware::buffers::SU3 * out, const hardware::buffers::Plain<Matrixsu3> * in);
	void convertGaugefieldFromSOA_device(const hardware::buffers::Plain<Matrixsu3> * out, const hardware::buffers::SU3 * in);

	/**
	 * Collect the kernels for OpenCL.
	 */
	void fill_kernels();
	/**
	 * Clear out the kernels,
	 */
	void clear_kernels();
};

#endif /* _OPENCL_MODULE_GAUGEFIELD_H_ */
