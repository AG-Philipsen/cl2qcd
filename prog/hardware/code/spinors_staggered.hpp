/** @file
 * Spinors_staggered class declaration
 */
#ifndef _HARDWARE_CODE_SPINORS_STAGGERED
#define _HARDWARE_CODE_SPINORS_STAGGERED

#include "opencl_module.hpp"

#include "../buffers/plain.hpp"
#include "../buffers/su3vec.hpp"
#include "../buffers/prng_buffer.hpp"

namespace hardware {

namespace code {

/**
 * An OpenCL device
 *
 * @todo Everything is public to faciliate inheritance. Actually, more parts should be private.
 */
class Spinors_staggered : public Opencl_Module {
public:
	friend hardware::Device;
	friend hardware::buffers::SU3vec;

	virtual ~Spinors_staggered();

	/*********************************************************************************************/
	/**************************  NON EVEN-ODD PRECONDITIONING METHODS  ***************************/
	/*********************************************************************************************/
	
	/////////////////////////////////////////////
	//        Fields algebra operations        //
	/////////////////////////////////////////////
	/**
	 * This function returns the input staggered field
	 * multiplied by a complex constant alpha
	 * @param x The input staggered field (one su3vec per site)
	 * @param alpha The complex constant
	 * @param out The output staggered field: alpha*x (one su3vec per site)
	 */
	void sax_device(const hardware::buffers::Plain<su3vec> * x, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Plain<su3vec> * out) const;
	
	/**
	 * This function returns the input staggered field x
	 * multiplied by alpha and added to the other input staggered field y
	 * @param x The first input staggered field (one su3vec per site)
	 * @param y The second input staggered field (one su3vec per site)
	 * @param alpha The complex constant
	 * @param out The output staggered field alpha*x + y (one su3vec per site)
	 */
	void saxpy_device(const hardware::buffers::Plain<su3vec> * x, const hardware::buffers::Plain<su3vec> * y, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Plain<su3vec> * out) const;
	
	/**
	 * This function returns a linear combination of the input staggered field
	 * @param x The first input staggered field (one su3vec per site)
	 * @param y The second input staggered field (one su3vec per site)
	 * @param z The third input staggered field (one su3vec per site)
	 * @param alpha The first complex constant
	 * @param beta The second complex constant
	 * @param out The output staggered field alpha*x + beta*y + z (one su3vec per site)
	 */
	void saxpbypz_device(const hardware::buffers::Plain<su3vec> * x, const hardware::buffers::Plain<su3vec> * y, const hardware::buffers::Plain<su3vec> * z, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Plain<hmc_complex> * beta, const hardware::buffers::Plain<su3vec> * out) const;
	
	//////////////////////////////////
	//      Algebra operations      //
	//////////////////////////////////
	/**
	 * This function completes the reduction to calculate the sum
	 * of squarenorms of the staggered field. It should be called
	 * after that the kernel "global_squarenorm_staggered" ended.
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
	
	/**
	 * This function initialize a staggered field with gaussian complex numbers
	 * @param in The field to be initialized 
	 * @param prng The cl_memory object
	 * @note A complex gaussian number is a number with both real and imaginary
	 *       part distributed in a gaussian way. In order to make the complex
	 *       gaussian distribution be normal (mean zero and variance 1), the real
	 *       and the imaginary parts must have mean zero and variance 0.5
	 *       It is worth remarking that a COMPLEX gaussian distribution is 
	 *       p(z)=1/(pi*sigma_z^2) exp(-|z|^2/sigma_z^2)
	 *       so it is a product of two real gaussian distributions on condition that
	 *       they have the SAME VARIANCE.
	 *       The complex variance sigma_z^2 will be both equal to 2*sigma_x^2 and
	 *       equal to 2*sigma_y. Hence to have sigma_z^2=1 we have to make
	 *       sigma_x^2=sigma_y^2=0.5 (for a very detailed treatment of this aspect
	 *       see "Statistical Signal Processing of Complex-Valued Data" from 
	 *       PETER J. SCHREIER around page 22).
	 */
	void set_gaussian_spinorfield_device(const hardware::buffers::Plain<su3vec> * in, const hardware::buffers::PRNGBuffer * prng) const;
	
	
	/*********************************************************************************************/
	/****************************  EVEN-ODD PRECONDITIONING METHODS  *****************************/
	/*********************************************************************************************/
	
	/////////////////////////////////////////////
	//        Fields algebra operations        //
	/////////////////////////////////////////////
	/**
	 * This function returns the input staggered field with even-odd preconditioning
	 * multiplied by a complex constant alpha
	 * @param x The input staggered field (one su3vec per site even or odd)
	 * @param alpha The complex constant
	 * @param out The output staggered field: alpha*x (one su3vec per site even or odd)
	 */
	void sax_eoprec_device(const hardware::buffers::SU3vec * x, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::SU3vec * out) const;
	
	/**
	 * This function returns the input staggered field x multiplied by alpha
	 * and added to the other input staggered field y (with even-odd preconditioning)
	 * @param x The first input staggered field (one su3vec per site even or odd)
	 * @param y The second input staggered field (one su3vec per site even or odd)
	 * @param alpha The complex constant
	 * @param out The output staggered field alpha*x + y (one su3vec per site even or odd)
	 */
	void saxpy_eoprec_device(const hardware::buffers::SU3vec * x, const hardware::buffers::SU3vec * y, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::SU3vec * out) const;	
	
	/**
	 * This function returns a linear combination of the input staggered fields
	 * (with even-odd preconditioning)
	 * @param x The first input staggered field (one su3vec per site even or odd)
	 * @param y The second input staggered field (one su3vec per site even or odd)
	 * @param alpha The first complex constant
	 * @param beta The second complex constant
	 * @param out The output staggered field alpha*x + beta*y (one su3vec per site even or odd)
	 */
	void saxpby_eoprec_device(const hardware::buffers::SU3vec * x, const hardware::buffers::SU3vec * y, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Plain<hmc_complex> * beta, const hardware::buffers::SU3vec * out) const;
	
	/**
	 * This function returns a linear combination of the input staggered fields
	 * (with even-odd preconditioning)
	 * @param x The first input staggered field (one su3vec per site even or odd)
	 * @param y The second input staggered field (one su3vec per site even or odd)
	 * @param z The third input staggered field (one su3vec per site even or odd)
	 * @param alpha The first complex constant
	 * @param beta The second complex constant
	 * @param out The output staggered field alpha*x + beta*y + z (one su3vec per site even or odd)
	 */
	void saxpbypz_eoprec_device(const hardware::buffers::SU3vec * x, const hardware::buffers::SU3vec * y, const hardware::buffers::SU3vec * z, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Plain<hmc_complex> * beta, const hardware::buffers::SU3vec * out) const;
	
	//////////////////////////////////
	//      Algebra operations      //
	//////////////////////////////////
	/**
	 * This function performs a complete reduction, calculating the sum of squarenorms
	 * of the staggered field with EVEN-ODD preconditioning.
	 * It uses the "global_squarenorm_reduction" function.
	 * @param a The staggered field (one su3vec per site)
	 * @param out The result of the reduction
	 */
	void set_float_to_global_squarenorm_eoprec_device(const hardware::buffers::SU3vec * a, const hardware::buffers::Plain<hmc_float> * out) const;
	
	/**
	 * This function calculates the complex scalarproduct of two staggered fields
	 * with even-odd preconditioning that means the sum of the complex scalar product
	 * of each su3vec per site (even or odd)
	 * @param a The first staggered field (one su3vec per site even or odd)
	 * @param b The second staggered field (one su3vec per site even or odd)
	 * @param out The result of the scalarproduct
	 */
	void set_complex_to_scalar_product_eoprec_device(const hardware::buffers::SU3vec * a, const hardware::buffers::SU3vec * b, const hardware::buffers::Plain<hmc_complex> * out) const;
	
	//////////////////////////////////
	//      Setting operations      //
	//////////////////////////////////
	/**
	 * This function sets to zero a staggered field with even-odd preconditioning.
	 * @param x The field to be set to zero 
	 */
	void set_zero_spinorfield_eoprec_device(const hardware::buffers::SU3vec * x) const;
	
	/**
	 * This function sets to 1 a staggered field and normalizes it (i.e. all its su3vec
	 * are set to 1./(sqrt(3*VOL4D_GLOBAL)). So the squarenorm of the field is 0.5 and not 1
	 * @param x The field to be set to one and normalized 
	 */
	void set_cold_spinorfield_eoprec_device(const hardware::buffers::SU3vec * x) const;
	
	/**
	 * This function initialize a staggered field with gaussian complex numbers
	 * (with even-odd preconditioning)
	 * @param in The field to be initialized 
	 * @param prng The cl_memory object
	 * @note A complex gaussian number is a number with both real and imaginary
	 *       part distributed in a gaussian way. In order to make the complex
	 *       gaussian distribution be normal (mean zero and variance 1), the real
	 *       and the imaginary parts must have mean zero and variance 0.5
	 *       It is worth remarking that a COMPLEX gaussian distribution is 
	 *       p(z)=1/(pi*sigma_z^2) exp(-|z|^2/sigma_z^2)
	 *       so it is a product of two real gaussian distributions on condition that
	 *       they have the SAME VARIANCE.
	 *       The complex variance sigma_z^2 will be both equal to 2*sigma_x^2 and
	 *       equal to 2*sigma_y. Hence to have sigma_z^2=1 we have to make
	 *       sigma_x^2=sigma_y^2=0.5 (for a very detailed treatment of this aspect
	 *       see "Statistical Signal Processing of Complex-Valued Data" from 
	 *       PETER J. SCHREIER around page 22).
	 */
	void set_gaussian_spinorfield_eoprec_device(const hardware::buffers::SU3vec * in, const hardware::buffers::PRNGBuffer * prng) const;
	
	/*********************************************************************************************/
	/************************************  GENERAL METHODS  **************************************/
	/*********************************************************************************************/
	
	//////////////////////////////////////////
	//     Conversions eo to/from non eo    //
	//////////////////////////////////////////
	/**
	 * This function reconstructs a field on the whole lattice from
	 * the fields on the even and on the odd sites
	 * @param even The staggered field on even sites
	 * @param odd The staggered field on odd sites
	 * @param out The staggered field on the whole lattice
	 */
	void convert_from_eoprec_device(const hardware::buffers::SU3vec * even, const hardware::buffers::SU3vec * odd, const hardware::buffers::Plain<su3vec> * out) const;
	
	/**
	 * This function splits a field on the whole lattice into
	 * two fields on the even and on the odd sites
	 * @param in The staggered field on the whole lattice to be splitted
	 * @param even The staggered field on even sites
	 * @param odd The staggered field on odd sites
	 */
	void convert_to_eoprec_device(const hardware::buffers::SU3vec * even, const hardware::buffers::SU3vec * odd, const hardware::buffers::Plain<su3vec> * in) const;
	
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

	//Functionalities to switch from AoS to SoA and viceversa
	cl_kernel convert_staggered_field_to_SoA_eo;
	cl_kernel convert_staggered_field_from_SoA_eo;

	/**
	 * This function reads a staggered field wrote in the AoS fashion
	 * and returns the same field wrote in the SoA fashion. Of course,
	 * in the SoA way, all first components of the su3vec are before the
	 * second, that are before the third.
	 * @param out The staggered field wrote in the SoA way
	 * @param in The staggered field wrote in the AoS way
	 * 
	 * @NOTE The input buffer is Plain and the output buffer will be SU3vec.
	 */
	void convert_staggered_field_to_SoA_eo_device(const hardware::buffers::SU3vec * out, const hardware::buffers::Plain<su3vec> * in) const;
	
	/**
	 * This function reads a staggered field wrote in the SoA fashion
	 * and returns the same field wrote in the AoS fashion. 
	 * @param out The staggered field wrote in the AoS way
	 * @param in The staggered field wrote in the SoA way
	 * 
	 * @NOTE The input buffer is of SU3vec type, while the output buffer will be Plain.
	 */
	void convert_staggered_field_from_SoA_eo_device(const hardware::buffers::Plain<su3vec> * out, const hardware::buffers::SU3vec * in) const;
	
	//Source for kernels
	ClSourcePackage basic_fermion_code;

	/******************************************************/
	/*******  NON EVEN-ODD PRECONDITIONING KERNELS  *******/
	/******************************************************/
	
	//Scalar Product
	cl_kernel scalar_product_stagg;
	cl_kernel scalar_product_reduction_stagg;
	
	//Squarenorm
	cl_kernel global_squarenorm_stagg;
	cl_kernel global_squarenorm_reduction_stagg;
	
	//Setting field
	cl_kernel set_zero_spinorfield_stagg;
	cl_kernel set_cold_spinorfield_stagg;
	cl_kernel set_gaussian_spinorfield_stagg;
	
	//Algebra on staggered fields
	cl_kernel sax_stagg;
	cl_kernel saxpy_stagg;
	cl_kernel saxpbypz_stagg;
	
	/******************************************************/
	/*******    EVEN-ODD PRECONDITIONING KERNELS    *******/
	/******************************************************/
	
	//Scalar Product
	cl_kernel scalar_product_stagg_eoprec;
	
	//Squarenorm
	cl_kernel global_squarenorm_stagg_eoprec;
	
	//Setting field
	cl_kernel set_zero_spinorfield_stagg_eoprec;
	cl_kernel set_cold_spinorfield_stagg_eoprec;
	cl_kernel set_gaussian_spinorfield_stagg_eoprec;
	
	//Algebra on staggered fields
	cl_kernel sax_stagg_eoprec;
	cl_kernel saxpy_stagg_eoprec;
	cl_kernel saxpby_stagg_eoprec;
	cl_kernel saxpbypz_stagg_eoprec;
	
	/******************************************************/
	/****************  GENERAL KERNELS  *******************/
	/******************************************************/
	
	//Conversion from even-odd preconditioning to non even-odd and viceversa
	cl_kernel convert_from_eoprec_stagg;
	cl_kernel convert_to_eoprec_stagg;
	
	
	/* To be added, if needed.
	
	//I would not implement these at first. So far, alpha as buffer everywhere.
	cl_kernel saxpy_arg_stagg; 
	cl_kernel saxpy_arg_stagg_eoprec;
	
	*/

};

}

}

#endif // _HARDWARE_CODE_SPINORS_
