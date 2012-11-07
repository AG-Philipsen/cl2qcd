/** @file
 * Basic OpenCL functionality
 */
#ifndef _OPENCLMODULEGAUGEMOMENTUMH_
#define _OPENCLMODULEGAUGEMOMENTUMH_

#include "hardware/code/opencl_module.hpp"
#include "hardware/buffers/plain.hpp"
#include "hardware/buffers/prng_buffer.hpp"
#include "hardware/buffers/gaugemomentum.hpp"

/**
 * An OpenCL device
 *
 * This class wraps all operations on a device. Operations are always specific, e.g. each kernel and copy operation
 * has it's own wrapper function.
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Opencl_Module_Gaugemomentum : public hardware::code::Opencl_Module {
public:
	friend hardware::Device;

	virtual ~Opencl_Module_Gaugemomentum();

	///////////////////////////////////////////////////
	//Methods on device
	void set_float_to_gaugemomentum_squarenorm_device(const hardware::buffers::Gaugemomentum * in, const hardware::buffers::Plain<hmc_float> * out);
	void generate_gaussian_gaugemomenta_device(const hardware::buffers::Gaugemomentum * in, const hardware::buffers::PRNGBuffer * prng);
	void set_zero_gaugemomentum(const hardware::buffers::Gaugemomentum *);
	/**
	 * Import data from the gaugemomenta array into the given buffer.
	 *
	 * The data in the buffer will be stored in the device specific format.
	 *
	 * @param[out] dest The buffer to write to in the device specific format
	 * @param[in] data The data to write to the buffer
	 */
	void importGaugemomentumBuffer(const hardware::buffers::Gaugemomentum * dest, const ae * const data);
	/**
	 * Export data from the given buffer into a normal gaugemomentum array.
	 *
	 * The data in the buffer is assumed to be in the device specific format.
	 *
	 * @param[out] dest An array that the buffer data can be written to.
	 * @param[in] data A buffer containing the data in the device specific format.
	 */
	void exportGaugemomentumBuffer(ae * const dest, const hardware::buffers::Gaugemomentum * buf);

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
	Opencl_Module_Gaugemomentum(const meta::Inputparameters& params, hardware::Device * device);

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

#endif //OPENCLMODULEGAUGEMOMENTUMH
