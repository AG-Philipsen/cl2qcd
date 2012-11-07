/** @file
 * Heatbath for OpenCL
 */
#ifndef _HARDWARE_CODE_SPINORS_
#define _HARDWARE_CODE_SPINORS_

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
class Spinors : public Opencl_Module {
public:
	friend hardware::Device;

	virtual ~Spinors();

	/////////////////////////////////////////
	// device operations
	void generate_gaussian_spinorfield_device(const hardware::buffers::Plain<spinor> * in, const hardware::buffers::PRNGBuffer * prng);
	void generate_gaussian_spinorfield_eo_device(const hardware::buffers::Spinor * in, const hardware::buffers::PRNGBuffer * prng);

	//    linear Algebra operations
	void convert_from_eoprec_device(const hardware::buffers::Spinor * in1, const hardware::buffers::Spinor * in2, const hardware::buffers::Plain<spinor> * out);
	void convert_to_eoprec_device(const hardware::buffers::Spinor * out1, const hardware::buffers::Spinor * out2, const hardware::buffers::Plain<spinor> * in);

	void set_complex_to_scalar_product_device(const hardware::buffers::Plain<spinor> * a, const hardware::buffers::Plain<spinor> * b, const hardware::buffers::Plain<hmc_complex> * out);
	void set_complex_to_scalar_product_eoprec_device(const hardware::buffers::Spinor * a, const hardware::buffers::Spinor * b, const hardware::buffers::Plain<hmc_complex> * out);
	void set_complex_to_ratio_device(const hardware::buffers::Plain<hmc_complex> * a, const hardware::buffers::Plain<hmc_complex> * b, const hardware::buffers::Plain<hmc_complex> * out);
	void set_complex_to_product_device(const hardware::buffers::Plain<hmc_complex> * a, const hardware::buffers::Plain<hmc_complex> * b, const hardware::buffers::Plain<hmc_complex> * out);
	void global_squarenorm_reduction(const hardware::buffers::Plain<hmc_float> * out, const hardware::buffers::Plain<hmc_float> * tmp_buf);
	void set_float_to_global_squarenorm_device(const hardware::buffers::Plain<spinor> * a, const hardware::buffers::Plain<hmc_float> * out);
	void set_float_to_global_squarenorm_eoprec_device(const hardware::buffers::Spinor * a, const hardware::buffers::Plain<hmc_float> * out);
	void set_zero_spinorfield_device(const hardware::buffers::Plain<spinor> * x);
	void set_zero_spinorfield_eoprec_device(const hardware::buffers::Spinor * x);
	void saxpy_device(const hardware::buffers::Plain<spinor> * x, const hardware::buffers::Plain<spinor> * y, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Plain<spinor> * out);
	void sax_device(const hardware::buffers::Plain<spinor> * x, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Plain<spinor> * out);
	void saxsbypz_device(const hardware::buffers::Plain<spinor> * x, const hardware::buffers::Plain<spinor> * y, const hardware::buffers::Plain<spinor> * z, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Plain<hmc_complex> * beta, const hardware::buffers::Plain<spinor> * out);
	void saxpy_eoprec_device(const hardware::buffers::Spinor * x, const hardware::buffers::Spinor * y, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Spinor * out);
	void sax_eoprec_device(const hardware::buffers::Spinor * x, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Spinor * out);
	void saxsbypz_eoprec_device(const hardware::buffers::Spinor * x, const hardware::buffers::Spinor * y, const hardware::buffers::Spinor * z, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Plain<hmc_complex> * beta, const hardware::buffers::Spinor * out);
	void create_point_source_device(const hardware::buffers::Plain<spinor> * inout, int i, int spacepos, int timepos);
	void create_point_source_eoprec_device(const hardware::buffers::Spinor * inout_even, const hardware::buffers::Spinor * inout_odd, cl_mem gf, int i, int spacepos, int timepos);
	void set_spinorfield_cold_device(const hardware::buffers::Plain<spinor> * inout);
	void set_eoprec_spinorfield_cold_device(const hardware::buffers::Spinor * inout);

	//    merged kernel calls
	void saxpy_AND_squarenorm_eo_device(const hardware::buffers::Spinor * x, const hardware::buffers::Spinor * y, const hardware::buffers::Plain<hmc_complex> * alpha, const hardware::buffers::Spinor * out, const hardware::buffers::Plain<hmc_complex> * sq_out);

	/**
	 * Copy an even-odd preconditioned spinorfield to the given buffer.
	 *
	 * @param buf A buffer of at least get_eoprec_spinorfield_buffer_size() bytes which will
	 *            be filled with the spinorfield in an implementation chosen format.
	 * @param source An array of spinors representing an even-odd field.
	 */
	void copy_to_eoprec_spinorfield_buffer(const hardware::buffers::Spinor * buf, const spinor * const source);

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
	Spinors(const meta::Inputparameters& params, hardware::Device * device);

	/**
	 * Collect the kernels for OpenCL.
	 */
	void fill_kernels();

	/**
	 * Clear out the kernels,
	 */
	void clear_kernels();

	cl_kernel convertSpinorfieldToSOA_eo;
	cl_kernel convertSpinorfieldFromSOA_eo;

	void convertSpinorfieldToSOA_eo_device(const hardware::buffers::Spinor * out, const hardware::buffers::Plain<spinor> * in);
	void convertSpinorfieldFromSOA_eo_device(const hardware::buffers::Plain<spinor> * out, const hardware::buffers::Spinor * in);

	ClSourcePackage basic_fermion_code;

	//BLAS
	cl_kernel set_spinorfield_cold;
	cl_kernel saxpy;
	cl_kernel sax;
	cl_kernel saxsbypz;
	cl_kernel set_zero_spinorfield;

	cl_kernel convert_from_eoprec;
	cl_kernel convert_to_eoprec;
	cl_kernel set_eoprec_spinorfield_cold;
	cl_kernel set_zero_spinorfield_eoprec;
	cl_kernel saxpy_eoprec;
	cl_kernel sax_eoprec;
	cl_kernel saxsbypz_eoprec;

	//Scalar Product
	cl_kernel scalar_product;
	cl_kernel scalar_product_reduction;
	cl_kernel global_squarenorm;
	cl_kernel _global_squarenorm_reduction;
	cl_kernel scalar_product_eoprec;
	cl_kernel global_squarenorm_eoprec;

	cl_kernel generate_gaussian_spinorfield;
	cl_kernel generate_gaussian_spinorfield_eo;

	//Single
	cl_kernel ratio;
	cl_kernel product;

	//merged kernels
	cl_kernel saxpy_AND_squarenorm_eo;
};

}

}

#endif // _HARDWARE_CODE_SPINORS_
