/** @file
 * Heatbath for OpenCL
 */
#ifndef _HARDWARE_CODE_SPINORS_STAGGERED
#define _HARDWARE_CODE_SPINORS_STAGGERED

#include "opencl_module.hpp"

#include "../buffers/plain.hpp"
#include "../buffers/spinor.hpp"
#include "../buffers/prng_buffer.hpp"

namespace hardware {

namespace code {

/**
 * An OpenCL device
 *
 * Adds random numbers to basic Opencl_Module class
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Spinors_staggered : public Opencl_Module {
public:
	friend hardware::Device;

	virtual ~Spinors_staggered();

	/////////////////////////////////////////
	//    linear Algebra operations        //
	/////////////////////////////////////////
	/**
	 * This function completes the reduction to calculate the sum
	 * of squarenorms of the staggered field. It should be called
	 * after that the kernel "global_squarenorm_staggered" ended.
	 * 
	 * @param out The result of the reduction
	 * @param tmp_buf The vector containing the results of the local (i.e. within each group)
	 *                reductions.
	 */
	void global_squarenorm_reduction(const hardware::buffers::Plain<hmc_float> * out, const hardware::buffers::Plain<hmc_float> * tmp_buf) const;
	/**
	 * This function performs a complete reduction, calculating the sum of squarenorms
	 * of the staggered field. It uses the "global_squarenorm_reduction" function.
	 * @param a The staggered field (one su3vec per site)
	 * @param out The result of the reduction
	 */
	void set_float_to_global_squarenorm_device(const hardware::buffers::Plain<su3vec> * a, const hardware::buffers::Plain<hmc_float> * out) const;

	///////////////////////////////////////////
	/**
	 * Print the profiling information to a file.
	 *
	 * @param filename Name of file where data is appended.
	 */
	void virtual print_profiling(const std::string& filename, int number) const override;

	ClSourcePackage get_sources() const noexcept;

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

protected:
	/**
	 * Add specific work_size determination for this child class
	 */
	virtual void get_work_sizes(const cl_kernel kernel, size_t * ls, size_t * gs, cl_uint * num_groups) const override;

private:
	Spinors_staggered(const meta::Inputparameters& params, hardware::Device * device);

	/**
	 * Collect the kernels for OpenCL.
	 */
	void fill_kernels();

	/**
	 * Clear out the kernels,
	 */
	void clear_kernels();

	ClSourcePackage basic_fermion_code;

	//Scalar Product
	cl_kernel global_squarenorm;
	cl_kernel _global_squarenorm_reduction;

};

}

}

#endif // _HARDWARE_CODE_SPINORS_
