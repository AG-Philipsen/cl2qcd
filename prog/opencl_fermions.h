/** @file
 * OpenCL device managment and everything executed on them -- including fermions.
 */
#ifndef _MYOPENCLFERMIONSH_
#define _MYOPENCLFERMIONSH_

#include "opencl.h"
//CP: this includes the struct-definitions for the spinors...
#include "types_fermions.h"
/** @todo is this needed?!?!? */
//#include "host_operations_spinorfield.h"

/**
 * An OpenCL device for fermionic calculations.
 *
 * This class wraps all operations on a device. Include fermions. Inherited from class Opencl.
 */
class Opencl_fermions : public Opencl {
public:
	/**
	* Collect a vector of kernel file names.
	* Virtual method, allows to include more kernel files in inherited classes.
	*/
	virtual hmc_error fill_kernels_file ();
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
	virtual hmc_error fill_kernels(cl_program program);


	/**
	 * Initialize the OpenCL device including fermion capabilities
	 *
	 * @param wanted The OpenCL device type to be used, e.g. CL_DEVICE_TYPE_CPU or CL_DEVICE_TYPE_GPU
	 * @param timer The timer to use for reporting execution time
	 * @param parameters The parsed input parameters
	 * @return Error code as defined in hmcerrs.h:
	 *         @li HMC_OCLERROR if OpenCL initialization / operations fail
	 *         @li HMC_FILEERROR if one of the kernel files cannot be opened
	 *         @li HMC_SUCCESS otherwise
	 */
	virtual hmc_error init(cl_device_type wanted_device_type, usetimer* timer, inputparameters* parameters);

	hmc_error finalize_fermions();

	////////////////////////////////////////////////7777
	// copying
	// now with the new spinor types
	//    non-eoprec
	hmc_error init_fermion_variables(inputparameters* parameters, usetimer* timer);
	hmc_error copy_spinorfield_to_device(spinorfield* host_spinorfield, usetimer* timer);
	hmc_error copy_source_to_device(spinorfield* host_spinorfield, usetimer* timer);
	hmc_error get_spinorfield_from_device(spinorfield* host_spinorfield,  usetimer* timer);
	hmc_error copy_spinor_device(cl_mem in, cl_mem out, usetimer* timer);

	//    eoprec
	hmc_error copy_eoprec_spinorfield_to_device(spinorfield_eoprec* host_spinorfield, usetimer* timer);
	hmc_error copy_eoprec_source_to_device(spinorfield_eoprec* host_spinorfield1, spinorfield_eoprec* host_spinorfield2, usetimer* timer);
	hmc_error get_eoprec_spinorfield_from_device(spinorfield_eoprec* host_spinorfield,  usetimer* timer);
	hmc_error copy_eoprec_spinor_device(cl_mem in, cl_mem out, usetimer* timer);

	//    misc
	hmc_error copy_float_from_device(cl_mem in, hmc_float * out, usetimer* timer);
	hmc_error copy_complex_from_device(cl_mem in, hmc_complex * out, usetimer* timer);
	hmc_error copy_complex_device(cl_mem in, cl_mem out, usetimer* timer);

	/////////////////////////////////////////
	// device operations

	//    linear Algebra operations
	hmc_error convert_to_kappa_format_device(cl_mem inout, const size_t local_work_size, const size_t global_work_size, usetimer* timer);
	hmc_error convert_from_kappa_format_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer);
	hmc_error convert_to_kappa_format_eoprec_device(cl_mem inout, const size_t local_work_size, const size_t global_work_size, usetimer* timer);
	hmc_error convert_from_kappa_format_eoprec_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer);
	hmc_error convert_from_eoprec_device(cl_mem in1, cl_mem in2, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer);

	hmc_error set_complex_to_scalar_product_device(cl_mem a, cl_mem b, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer);
	hmc_error set_complex_to_scalar_product_eoprec_device(cl_mem a, cl_mem b, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer);
	hmc_error set_complex_to_ratio_device(cl_mem a, cl_mem b, cl_mem out, usetimer* timer);
	hmc_error set_complex_to_product_device(cl_mem a, cl_mem b, cl_mem out, usetimer* timer);
	hmc_error set_float_to_global_squarenorm_device(cl_mem a, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer);
	hmc_error set_float_to_global_squarenorm_eoprec_device(cl_mem a, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer);
	hmc_error set_zero_spinorfield_device(cl_mem x, const size_t local_work_size, const size_t global_work_size, usetimer* timer);
	hmc_error set_zero_spinorfield_eoprec_device(cl_mem x, const size_t local_work_size, const size_t global_work_size, usetimer* timer);
	hmc_error saxpy_device(cl_mem x, cl_mem y, cl_mem alpha, cl_mem out, const size_t local_work_size, const size_t global_work_size,  usetimer* timer);
	hmc_error saxsbypz_device(cl_mem x, cl_mem y, cl_mem z, cl_mem alpha, cl_mem beta, cl_mem out, const size_t local_work_size, const size_t global_work_size,  usetimer* timer);
	hmc_error saxpy_eoprec_device(cl_mem x, cl_mem y, cl_mem alpha, cl_mem out, const size_t local_work_size, const size_t global_work_size,  usetimer* timer);
	hmc_error saxsbypz_eoprec_device(cl_mem x, cl_mem y, cl_mem z, cl_mem alpha, cl_mem beta, cl_mem out, const size_t local_work_size, const size_t global_work_size,  usetimer* timer);
	hmc_error create_point_source_device(int i, int spacepos, int timepos, const size_t ls, const size_t gs, usetimer * latimer);
	hmc_error create_point_source_eoprec_device(int i, int spacepos, int timepos, const size_t ls, const size_t gs, usetimer * latimer, usetimer * dslashtimer, usetimer * Mdiagtimer);
	hmc_error set_spinorfield_cold_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer);
	hmc_error set_eoprec_spinorfield_cold_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer);

	//    fermionmatrix operations
	//    non-eoprec
	hmc_error M_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size,  usetimer* timer, usetimer * dslashtimer, usetimer * Mdiagtimer);
	hmc_error gamma5_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer *timer);
	hmc_error Qplus_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size,  usetimer* timer);
	hmc_error Qminus_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size,  usetimer* timer);
	hmc_error QplusQminus_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size,  usetimer* timer);
	//    eoprec
	hmc_error gamma5_eoprec_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer *timer);
	hmc_error Aee_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer, usetimer * singletimer, usetimer * dslashtimer, usetimer * Mdiagtimer, usetimer * latimer);
	hmc_error M_inverse_sitediagonal_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer * timer);
	hmc_error M_sitediagonal_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer * timer);
	hmc_error dslash_eoprec_device(cl_mem in, cl_mem out, int evenodd, const size_t local_work_size, const size_t global_work_size, usetimer * timer);

	//this is not needed anymore!!
//  hmc_error testing_spinor(inputparameters* parameters, size_t local_size, size_t global_size);

	//    solver operations
	//    non-eoprec
	hmc_error solver_device(usetimer * copytimer, usetimer * singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, usetimer * dslashtimer, usetimer * Mdiagtimer, usetimer * solvertimer, const size_t ls, const size_t gs, int cgmax);
	hmc_error bicgstab_device(usetimer * copytimer, usetimer* singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, usetimer * dslashtimer, usetimer * Mdiagtimer, const size_t local_work_size, const size_t global_work_size, int cgmax);
	hmc_error cg_device(usetimer * copytimer, usetimer* singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, usetimer * dslashtimer, usetimer * Mdiagtimer, const size_t local_work_size, const size_t global_work_size, int cgmax);
	//    eorec
	hmc_error bicgstab_eoprec_device(usetimer * copytimer, usetimer* singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, usetimer * dslashtimer, usetimer * Mdiagtimer, const size_t local_work_size, const size_t global_work_size, int cgmax);
	hmc_error solver_eoprec_device(usetimer * copytimer, usetimer * singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, usetimer * dslashtimer, usetimer * Mdiagtimer, usetimer * solvertimer, const size_t ls, const size_t gs, int cgmax);

	//    operations needed calculating fermionic observables
	hmc_error ps_correlator_device(const size_t local_work_size, const size_t global_work_size, usetimer * timer);
	hmc_error set_correlator_field_zero_device(const size_t ls, const size_t gs, usetimer * latimer);
	hmc_error add_solution_to_correlator_field_device(const size_t ls, const size_t gs, usetimer * latimer);

	/////////////////////////////////////////////////
	//functions to get private variables
	cl_mem get_clmem_inout();
	cl_mem get_clmem_tmp();

	cl_mem get_clmem_inout_eoprec();
	cl_mem get_clmem_tmp_eoprec_1();

private:
	////////////////////////////////////
	//kernels, sorted roughly by groups
	//fermionmatrix
	cl_kernel M;
	cl_kernel gamma5;
	cl_kernel Qplus;
	cl_kernel Qminus;
	cl_kernel gamma5_eoprec;
	cl_kernel M_sitediagonal;
	cl_kernel M_inverse_sitediagonal;
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
