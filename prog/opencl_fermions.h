/** @file
 * OpenCL device managment and everything executed on them -- including fermions.
 */
#ifndef _MYOPENCLFERMIONSH_
#define _MYOPENCLFERMIONSH_


class Opencl_fermions;

#include "opencl.h"
//CP: this includes the struct-definitions for the spinors...
#include "types_fermions.h"
/** @todo is this needed?!?!? */
//#include "host_operations_spinorfield.h"

/**
 * this is a workaround to be able to pass a (fermionmatrix-)function, which are methods in this class,
 * to another function inside this class.
 * This type points to a helper-function, which then calls the wanted function.
 */
typedef hmc_error (*matrix_function_call) (Opencl_fermions* that, cl_mem in, cl_mem out, cl_mem gf, const size_t ls, const size_t gs);
hmc_error M_call(Opencl_fermions* that, cl_mem in, cl_mem out, cl_mem gf, const size_t ls, const size_t gs);
hmc_error Qplus_call(Opencl_fermions* that, cl_mem in, cl_mem out, cl_mem gf, const size_t ls, const size_t gs);
hmc_error Qminus_call(Opencl_fermions* that, cl_mem in, cl_mem out, cl_mem gf, const size_t ls, const size_t gs);
hmc_error QplusQminus_call(Opencl_fermions* that, cl_mem in, cl_mem out, cl_mem gf, const size_t ls, const size_t gs);
hmc_error Aee_call(Opencl_fermions* that, cl_mem in, cl_mem out, cl_mem gf, const size_t ls, const size_t gs);

/**
 * An OpenCL device for fermionic calculations.
 *
 * This class wraps all operations on a device. Include fermions. Inherited from class Opencl.
 */
class Opencl_fermions : public Opencl {
public:
	/**
	 * Collect the compiler options for OpenCL.
	 * Virtual method, allows to include more options in inherited classes.
	 */
	virtual hmc_error fill_collect_options(stringstream* collect_options);
	/**
	 * Collect the buffers to generate for OpenCL.
	 * Virtual method, allows to include more buffers in inherited classes.
	 */
	virtual hmc_error fill_buffers();
	/**
	 * Collect the kernels for OpenCL.
	 * Virtual method, allows to include more kernels in inherited classes.
	 */
	virtual void fill_kernels();

	/**
	 * Initialize the OpenCL device including fermion capabilities
	 *
	 * @param wanted The OpenCL device type to be used, e.g. CL_DEVICE_TYPE_CPU or CL_DEVICE_TYPE_GPU
	 * @param parameters The parsed input parameters
	 * @return Error code as defined in hmcerrs.h:
	 *         @li HMC_OCLERROR if OpenCL initialization / operations fail
	 *         @li HMC_FILEERROR if one of the kernel files cannot be opened
	 *         @li HMC_SUCCESS otherwise
	 */
	virtual hmc_error init(cl_device_type wanted_device_type, inputparameters* parameters);

	hmc_error finalize_fermions();

	////////////////////////////////////////////////
	hmc_error init_fermion_variables(inputparameters* parameters, usetimer* timer);

	/////////////////////////////////////////
	// device operations

	//    linear Algebra operations
	hmc_error convert_to_kappa_format_device(cl_mem inout, const size_t ls, const size_t gs);
	hmc_error convert_from_kappa_format_device(cl_mem in, cl_mem out, const size_t ls, const size_t gs);
	hmc_error convert_to_kappa_format_eoprec_device(cl_mem inout, const size_t ls, const size_t gs);
	hmc_error convert_from_kappa_format_eoprec_device(cl_mem in, cl_mem out, const size_t ls, const size_t gs);
	hmc_error convert_from_eoprec_device(cl_mem in1, cl_mem in2, cl_mem out, const size_t ls, const size_t gs);

	hmc_error set_complex_to_scalar_product_device(cl_mem a, cl_mem b, cl_mem out, const size_t ls, const size_t gs);
	hmc_error set_complex_to_scalar_product_eoprec_device(cl_mem a, cl_mem b, cl_mem out, const size_t ls, const size_t gs);
	hmc_error set_complex_to_ratio_device(cl_mem a, cl_mem b, cl_mem out);
	hmc_error set_complex_to_product_device(cl_mem a, cl_mem b, cl_mem out);
	hmc_error set_float_to_global_squarenorm_device(cl_mem a, cl_mem out, const size_t ls, const size_t gs);
	hmc_error set_float_to_global_squarenorm_eoprec_device(cl_mem a, cl_mem out, const size_t ls, const size_t gs);
	hmc_error set_zero_spinorfield_device(cl_mem x, const size_t ls, const size_t gs);
	hmc_error set_zero_spinorfield_eoprec_device(cl_mem x, const size_t ls, const size_t gs);
	hmc_error saxpy_device(cl_mem x, cl_mem y, cl_mem alpha, cl_mem out, const size_t ls, const size_t gs);
	hmc_error saxsbypz_device(cl_mem x, cl_mem y, cl_mem z, cl_mem alpha, cl_mem beta, cl_mem out, const size_t ls, const size_t gs);
	hmc_error saxpy_eoprec_device(cl_mem x, cl_mem y, cl_mem alpha, cl_mem out, const size_t ls, const size_t gs);
	hmc_error saxsbypz_eoprec_device(cl_mem x, cl_mem y, cl_mem z, cl_mem alpha, cl_mem beta, cl_mem out, const size_t ls, const size_t gs);
	void create_point_source_device(cl_mem inout, int i, int spacepos, int timepos, const size_t ls, const size_t gs);
	void create_point_source_eoprec_device(cl_mem inout_even, cl_mem inout_odd, cl_mem gf, int i, int spacepos, int timepos, const size_t ls, const size_t gs);
	void set_spinorfield_cold_device(cl_mem inout, const size_t ls, const size_t gs);
	void set_eoprec_spinorfield_cold_device(cl_mem inout, const size_t ls, const size_t gs);

	//    fermionmatrix operations
	//    non-eoprec
	//        compound
	hmc_error M(cl_mem in, cl_mem out, cl_mem gf, const size_t ls, const size_t gs);
	hmc_error Qplus(cl_mem in, cl_mem out, cl_mem gf, const size_t ls, const size_t gs);
	hmc_error Qminus(cl_mem in, cl_mem out, cl_mem gf, const size_t ls, const size_t gs);
	hmc_error QplusQminus(cl_mem in, cl_mem out, cl_mem gf, const size_t ls, const size_t gs);
	//        explicit
	hmc_error M_wilson_device(cl_mem in, cl_mem out, cl_mem gf, const size_t ls, const size_t gs);
	hmc_error M_tm_plus_device(cl_mem in, cl_mem out, cl_mem gf, const size_t ls, const size_t gs);
	hmc_error M_tm_minus_device(cl_mem in, cl_mem out, cl_mem gf, const size_t ls, const size_t gs);
	hmc_error gamma5_device(cl_mem inout, const size_t ls, const size_t gs);
	//    eoprec
	//        compound
	hmc_error Aee(cl_mem in, cl_mem out, cl_mem gf, const size_t ls, const size_t gs);
	//        explicit
	hmc_error gamma5_eoprec_device(cl_mem inout, const size_t ls, const size_t gs);
	hmc_error M_tm_inverse_sitediagonal_device(cl_mem in, cl_mem out, const size_t ls, const size_t gs);
	hmc_error M_tm_sitediagonal_device(cl_mem in, cl_mem out, const size_t ls, const size_t gs);
	hmc_error dslash_eoprec_device(cl_mem in, cl_mem out, cl_mem gf, int evenodd, const size_t ls, const size_t gs);

	//    solver operations
	//    non-eoprec
	/// this calls the solver according to parameter settings using the fermionmatrix f
	hmc_error solver_device(cl_mem gf, usetimer * copytimer, usetimer * singletimer, usetimer * solvertimer, const size_t ls, const size_t gs, int cgmax, matrix_function_call f);
	/// this executes the bicgstab on the device, using the fermionmatrix f
	hmc_error bicgstab_device(cl_mem gf, usetimer * copytimer, usetimer* singletimer, const size_t ls, const size_t gs, int cgmax, matrix_function_call f);
	/// this executes the cg on the device, using the fermionmatrix f 
	hmc_error cg_device(cl_mem gf, usetimer * copytimer, usetimer* singletimer, const size_t ls, const size_t gs, int cgmax, matrix_function_call f);
	//    eoprec
	/// this calls the solver according to parameter settings using the fermionmatrix f
	hmc_error solver_eoprec_device(cl_mem gf, usetimer * copytimer, usetimer * singletimer, usetimer * solvertimer, const size_t ls, const size_t gs, int cgmax, matrix_function_call f);	
	/// this executes the eoprec bicgstab on the device, using the fermionmatrix f 
	hmc_error bicgstab_eoprec_device(cl_mem gf, usetimer * copytimer, usetimer* singletimer, const size_t ls, const size_t gs, int cgmax, matrix_function_call f);
	hmc_error cg_eoprec_device(cl_mem gf, usetimer * copytimer, usetimer* singletimer, const size_t ls, const size_t gs, int cgmax, matrix_function_call f);
	
	//    operations needed calculating fermionic observables
	hmc_error ps_correlator_device(const size_t ls, const size_t gs);
	hmc_error set_correlator_field_zero_device(const size_t ls, const size_t gs);
	hmc_error add_solution_to_correlator_field_device(const size_t ls, const size_t gs);

	/////////////////////////////////////////////////
	//functions to get private variables
	cl_mem get_clmem_inout();
	cl_mem get_clmem_source();
	cl_mem get_clmem_tmp();
	cl_mem get_clmem_corr();

	cl_mem get_clmem_inout_eoprec();
	cl_mem get_clmem_tmp_eoprec_1();
	cl_mem get_clmem_source_even();
	cl_mem get_clmem_source_odd();

#ifdef _PROFILING_
	//CP: if PROFILING is activated, one needs a timer for each kernel
	//fermionmatrix
	usetimer timer_M_wilson;
	usetimer timer_gamma5;
	usetimer timer_M_tm_plus;
	usetimer timer_M_tm_minus;
	usetimer timer_gamma5_eoprec;
	usetimer timer_M_tm_sitediagonal;
	usetimer timer_M_tm_inverse_sitediagonal;
	usetimer timer_dslash_eoprec;

	//BLAS
	usetimer timer_set_spinorfield_cold;
	usetimer timer_set_eoprec_spinorfield_cold;
	usetimer timer_convert_from_eoprec;
	usetimer timer_saxpy;
	usetimer timer_saxsbypz;
	usetimer timer_set_zero_spinorfield;
	usetimer timer_convert_to_kappa_format;
	usetimer timer_convert_from_kappa_format;
	usetimer timer_convert_to_kappa_format_eoprec;
	usetimer timer_convert_from_kappa_format_eoprec;
	usetimer timer_create_point_source;
	usetimer timer_saxpy_eoprec;
	usetimer timer_saxsbypz_eoprec;
	usetimer timer_set_zero_spinorfield_eoprec;
	usetimer timer_create_point_source_eoprec;

	//Scalar Product
	usetimer timer_scalar_product;
	usetimer timer_scalar_product_reduction;
	usetimer timer_global_squarenorm;
	usetimer timer_global_squarenorm_reduction;
	usetimer timer_scalar_product_eoprec;
	usetimer timer_global_squarenorm_eoprec;

	//Single
	usetimer timer_ratio;
	usetimer timer_product;

	//Observables
	usetimer timer_ps_correlator;

	/**
	 * Return the timer connected to a specific kernel.
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual usetimer* get_timer(char * in);
	
	/**
	 * Return amount of bytes read and written by a specific kernel per call. 
	 *
	 * @param in Name of the kernel under consideration.
	 */
	virtual int get_read_write_size(char * in, inputparameters * parameters);	
	
	/**
	 * Print the profiling information to a file.
	 *
	 * @param filename Name of file where data is appended.
	 * @param parameters inputparameters
	 */
	void virtual print_profiling(std::string filename);	
#endif	
	
private:
	////////////////////////////////////
	//kernels, sorted roughly by groups
	//fermionmatrix
	cl_kernel M_wilson;
	cl_kernel gamma5;
	cl_kernel M_tm_plus;
	cl_kernel M_tm_minus;
	cl_kernel gamma5_eoprec;
	cl_kernel M_tm_sitediagonal;
	cl_kernel M_tm_inverse_sitediagonal;
	cl_kernel dslash_eoprec;

	//BLAS
	cl_kernel set_spinorfield_cold;
	cl_kernel set_eoprec_spinorfield_cold;
	cl_kernel convert_from_eoprec;
	cl_kernel saxpy;
	cl_kernel saxsbypz;
	cl_kernel set_zero_spinorfield;
	cl_kernel convert_to_kappa_format;
	cl_kernel convert_from_kappa_format;
	cl_kernel convert_to_kappa_format_eoprec;
	cl_kernel convert_from_kappa_format_eoprec;
	cl_kernel create_point_source;
	cl_kernel saxpy_eoprec;
	cl_kernel saxsbypz_eoprec;
	cl_kernel set_zero_spinorfield_eoprec;
	cl_kernel create_point_source_eoprec;

	//Scalar Product
	cl_kernel scalar_product;
	cl_kernel scalar_product_reduction;
	cl_kernel global_squarenorm;
	cl_kernel global_squarenorm_reduction;
	cl_kernel scalar_product_eoprec;
	cl_kernel global_squarenorm_eoprec;

	//Single
	cl_kernel ratio;
	cl_kernel product;

	//Observables
	cl_kernel ps_correlator;

	//CP: variables for normal solver
	cl_mem clmem_inout;
	cl_mem clmem_source;
	cl_mem clmem_rn;
	cl_mem clmem_rhat;
	cl_mem clmem_v;
	cl_mem clmem_p;
	cl_mem clmem_s;
	cl_mem clmem_t;
	cl_mem clmem_aux;
	//this is needed in QplusQminus as a temporary field
	cl_mem clmem_tmp;

	//CP: variables for eoprec solver
	cl_mem clmem_inout_eoprec;
	cl_mem clmem_source_even;
	cl_mem clmem_source_odd;
	cl_mem clmem_rn_eoprec;
	cl_mem clmem_rhat_eoprec;
	cl_mem clmem_v_eoprec;
	cl_mem clmem_p_eoprec;
	cl_mem clmem_s_eoprec;
	cl_mem clmem_t_eoprec;
	cl_mem clmem_aux_eoprec;
	cl_mem clmem_tmp_eoprec_1;
	cl_mem clmem_tmp_eoprec_2;
	cl_mem clmem_tmp_eoprec_3;

	cl_mem clmem_rho;
	cl_mem clmem_rho_next;
	cl_mem clmem_alpha;
	cl_mem clmem_omega;
	cl_mem clmem_beta;
	cl_mem clmem_tmp1;
	cl_mem clmem_tmp2;
	cl_mem clmem_one;
	cl_mem clmem_minusone;

	cl_mem clmem_kappa_cmplx;// = {kappa, 0.};
	cl_mem clmem_scalar_product_buf_glob;
	cl_mem clmem_global_squarenorm_buf_glob;
	cl_mem clmem_resid;
	cl_mem clmem_trueresid;

	//CP: variables for correlator, these should propably not be here if no observables should be calculated...
	cl_mem clmem_corr;

protected:
	ClSourcePackage basic_fermion_code;
};

#endif // _MYOPENCLFERMIONSH_
