/** @file
 * Basic OpenCL functionality
 */
#ifndef _HARDWARE_CODE_GAUGEMOMENTUM_
#define _HARDWARE_CODE_GAUGEMOMENTUM_

#include "opencl_module.hpp"
#include "../buffers/plain.hpp"
#include "../buffers/prng_buffer.hpp"
#include "../buffers/gaugemomentum.hpp"

namespace hardware {

namespace code {

/**
 * An OpenCL device
 *
 * This class wraps all operations on a device. Operations are always specific, e.g. each kernel and copy operation
 * has it's own wrapper function.
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Gaugemomentum : public Opencl_Module {
public:
	friend hardware::Device;

	virtual ~Gaugemomentum();

	///////////////////////////////////////////////////
	//Methods on device
	void set_float_to_gaugemomentum_squarenorm_device(const hardware::buffers::Gaugemomentum * in, const hardware::buffers::Plain<hmc_float> * out) const;
	void generate_gaussian_gaugemomenta_device(const hardware::buffers::Gaugemomentum * in, const hardware::buffers::PRNGBuffer * prng) const;
	void set_zero_gaugemomentum(const hardware::buffers::Gaugemomentum *) const;
	/**
	 * Import data from the gaugemomenta array into the given buffer.
	 *
	 * The data in the buffer will be stored in the device specific format.
	 *
	 * @param[out] dest The buffer to write to in the device specific format
	 * @param[in] data The data to write to the buffer
	 */
	void importGaugemomentumBuffer(const hardware::buffers::Gaugemomentum * dest, const ae * const data) const;
	/**
	 * Export data from the given buffer into a normal gaugemomentum array.
	 *
	 * The data in the buffer is assumed to be in the device specific format.
	 *
	 * @param[out] dest An array that the buffer data can be written to.
	 * @param[in] data A buffer containing the data in the device specific format.
	 */
	void exportGaugemomentumBuffer(ae * const dest, const hardware::buffers::Gaugemomentum * buf) const;

	ClSourcePackage get_sources() const noexcept;

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
	Gaugemomentum(const meta::Inputparameters& params, hardware::Device * device);

	/**
	 * Collect the kernels for OpenCL.
	 */
	void fill_kernels();
	/**
	 * Clear out the kernels,
	 */
	void clear_kernels();

	ClSourcePackage basic_gaugemomentum_code;

	//kernels
	cl_kernel generate_gaussian_gaugemomenta;
	cl_kernel _set_zero_gaugemomentum;
	cl_kernel gaugemomentum_squarenorm;
	cl_kernel gaugemomentum_convert_to_soa;
	cl_kernel gaugemomentum_convert_from_soa;
};

}

}

#endif // _HARDWARE_CODE_GAUGEMOMENTUM_
