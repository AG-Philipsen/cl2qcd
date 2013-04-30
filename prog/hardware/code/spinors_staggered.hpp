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

	/**
	 * This function calculates the complex scalarproduct of two staggered fields
	 * that means the sum of the complex scalar product of each su3vec per site
	 * @param a The first staggered field (one su3vec per site)
	 * @param b The second staggered field (one su3vec per site)
	 * @param out The result of the scalarproduct
	 */
	void set_complex_to_scalar_product_device(const hardware::buffers::Plain<su3vec> * a, const hardware::buffers::Plain<su3vec> * b, const hardware::buffers::Plain<hmc_complex> * out) const;
	
	//////////////////////////////////
	//      Setting operations      //
	//////////////////////////////////
	
	/**
	 * This function sets to zero a staggered field (all its su3vec)
	 * @param x The field to be set to zero 
	 */
	void set_zero_spinorfield_device(const hardware::buffers::Plain<su3vec> * x) const;
	
	/**
	 * This function sets to 1 a staggered field and normalizes it (i.e. all its su3vec
	 * are set to 1./(sqrt(3*VOL4D_GLOBAL))
	 * @param x The field to be set to one and normalized 
	 */
	void set_cold_spinorfield_device(const hardware::buffers::Plain<su3vec> * x) const;
	
	//////////////////////////////////////////
	//      Complex numbers operations      //
	//////////////////////////////////////////
	
	/**
	 * This function converts a float number into a complex one
	 * @param in The float number to be converted
	 * @param out The output complex number
	 */
	void set_complex_to_float_device(const hardware::buffers::Plain<hmc_float> * in, const hardware::buffers::Plain<hmc_complex> * out) const;
	
	/**
	 * This function executes the division between two complex numbers
	 * @param a The numerator
	 * @param b The denominator
	 * @param out Complex number a/b
	 */
	void set_complex_to_ratio_device(const hardware::buffers::Plain<hmc_complex> * a, const hardware::buffers::Plain<hmc_complex> * b, const hardware::buffers::Plain<hmc_complex> * out) const;
	
	/**
	 * This function executes the multiplication between two complex numbers
	 * @param a The first factor
	 * @param b The second factor
	 * @param out Complex number a*b
	 */
	void set_complex_to_product_device(const hardware::buffers::Plain<hmc_complex> * a, const hardware::buffers::Plain<hmc_complex> * b, const hardware::buffers::Plain<hmc_complex> * out) const;
	
	
	////////////////////////////////////////////////////////////////////////////////////////
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
	cl_kernel scalar_product_stagg;
	cl_kernel scalar_product_reduction_stagg;
	cl_kernel global_squarenorm_stagg;
	cl_kernel global_squarenorm_reduction_stagg;
	
	//Setting field
	cl_kernel set_zero_spinorfield_stagg;
	cl_kernel set_cold_spinorfield_stagg;
	
	//Operations between complex numbers
	cl_kernel convert_stagg;
	cl_kernel ratio_stagg;
	cl_kernel product_stagg;
	
	/* To be added...
	
	cl_kernel sax;
	cl_kernel saxpy;
	cl_kernel saxpy_arg;
	cl_kernel saxsbypz;
	cl_kernel generate_gaussian_spinorfield;
	
	*/

};

}

}

#endif // _HARDWARE_CODE_SPINORS_
