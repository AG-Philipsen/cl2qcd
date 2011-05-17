/** @file
 * OpenCL device managment and everything executed on them.
 */
#ifndef _MYOPENCLH_
#define _MYOPENCLH_

#include <cstdlib>
#include <vector>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include "host_geometry.h"
#include "host_operations_complex.h"
#include "host_operations_gaugefield.h"
#include "host_operations_spinor.h"
#include "host_operations_spinorfield.h"
#include "host_operations_fermionmatrix.h"
#include "globaldefs.h"
#include "hmcerrs.h"
#include "types.h"
#include "host_use_timer.h"
#include "host_testing.h"
#include "host_random.h"

/**
 * An OpenCL device
 *
 * This class wraps all operations on a device. Operations are always specific, e.g. each kernel and copy operation
 * has it's own wrapper function.
 */
class opencl {
public:
	/**
	 * Default constructor that also initializes the device.
	 *
	 * @param wanted The OpenCL device type to be used, e.g. CL_DEVICE_TYPE_CPU or CL_DEVICE_TYPE_GPU
	 * @param ls The local work size to be used on the device
	 * @param gs The global work size to be used on the device
	 * @param timer The timer to use for reporting execution time
	 * @param parameters The parsed input parameters
	 *
	 * @todo Should probably throw an exception on error
	 */
	opencl(cl_device_type wanted, const size_t ls, const size_t gs, usetimer* timer, inputparameters* parameters) {
		init(wanted, ls, gs, timer, parameters);
	};
	/**
	 * Empty constructor. Needed for gaugefield class.
	 */
	opencl() {};
	~opencl() {
		finalize();
	};

	/**
	 * Initialize the OpenCL device
	 *
	 * @param wanted The OpenCL device type to be used, e.g. CL_DEVICE_TYPE_CPU or CL_DEVICE_TYPE_GPU
	 * @param ls The local work size to be used on the device
	 * @param gs The global work size to be used on the device
	 * @param timer The timer to use for reporting execution time
	 * @param parameters The parsed input parameters
	 * @return Error code as defined in hmcerrs.h:
	 *         @li HMC_OCLERROR if OpenCL initialization / operations fail
	 *         @li HMC_FILEERROR if one of the kernel files cannot be opened
	 *         @li HMC_SUCCESS otherwise
	 */
	hmc_error init(cl_device_type wanted_device_type, const size_t local_work_size, const size_t global_work_size, usetimer* timer, inputparameters* parameters);

	/**
	 * Copy the given gaugefield to the appropriate OpenCL buffer.
	 *
	 * @param host_gaugefield The gaugefield to copy
	 * @param timer The timer to use to measure the copying time
	 * @return Error code as defined in hmcerrs.h:
	 *         @li HMC_OCLERROR if OpenCL operations fail
	 *         @li HMC_SUCCESS otherwise
	 */
	hmc_error copy_gaugefield_to_device(hmc_gaugefield* host_gaugefield,  usetimer* timer);

	/**
	 * Copy the RNG state to the appropriate OpenCL buffer.
	 *
	 * @param host_rndarray The RNG state to copy
	 * @param timer The timer to use to measure the copying time
	 * @return Error code as defined in hmcerrs.h:
	 *         @li HMC_OCLERROR if OpenCL operations fail
	 *         @li HMC_SUCCESS otherwise
	 */
	hmc_error copy_rndarray_to_device(hmc_rndarray host_rndarray,  usetimer* timer);

	hmc_error copy_rndarray_from_device(hmc_rndarray rndarray, usetimer* timer);

	/**
	 * Copy the gaugefield from the device into the given memory location.
	 *
	 * @param host_gaugefield Storage location for the gaugefield
	 * @param timer The timer to use to measure the copying time
	 * @return Error code as defined in hmcerrs.h:
	 *         @li HMC_OCLERROR if OpenCL operations fail
	 *         @li HMC_SUCCESS otherwise
	 */
	hmc_error get_gaugefield_from_device(hmc_gaugefield* host_gaugefield,  usetimer* timer);

	/**
	 * Perform one heatbath step.
	 */
	hmc_error run_heatbath(const hmc_float beta, usetimer * const timer);

	/**
	 * Perform one overrelaxation step.
	 */
	hmc_error run_overrelax(const hmc_float beta, usetimer * const timer);

	/**
	 * Calculate plaquette and polyakov.
	 *
	 * @param[in] local_work_size The local work size to use
	 * @param[in] global_work_size The global works size to use
	 * @param[out] plaq Storage for result of plaquette calculation
	 * @param[out] tplaq Storage for result of plaquette calculation
	 * @param[out] splaq Storage for result of plaquette calculation
	 * @param[out] pol Storage for result of polyakov calculation
	 * @param[in,out] timer1 Timer into which to aggregate plaquette calculation time
	 * @param[in,out] timer2 Timer into which to aggregate polyakov calculation time
	 * @return Error code as defined in hmcerrs.h:
	 *         @li HMC_OCLERROR if OpenCL operations fail
	 *         @li HMC_SUCCESS otherwise
	 */
	hmc_error gaugeobservables(const size_t local_work_size, const size_t global_work_size, hmc_float * plaq, hmc_float * tplaq, hmc_float * splaq, hmc_complex * pol, usetimer* timer1, usetimer* timer2);

#ifdef _FERMIONS_
	hmc_error init_fermion_variables(inputparameters* parameters, const size_t local_work_size, const size_t global_work_size, usetimer* timer);
	hmc_error copy_spinorfield_to_device(hmc_spinor_field* host_spinorfield, usetimer* timer);
	hmc_error copy_eoprec_spinorfield_to_device(hmc_spinor_field* host_spinorfield, usetimer* timer);
	hmc_error copy_source_to_device(hmc_spinor_field* host_spinorfield, usetimer* timer);
	hmc_error copy_eoprec_source_to_device(hmc_eoprec_spinor_field* host_spinorfield1, hmc_eoprec_spinor_field* host_spinorfield2, usetimer* timer);
	hmc_error get_spinorfield_from_device(hmc_spinor_field* host_spinorfield,  usetimer* timer);
	hmc_error get_eoprec_spinorfield_from_device(hmc_spinor_field* host_spinorfield,  usetimer* timer);
	hmc_error copy_spinor_device(cl_mem in, cl_mem out, usetimer* timer);
	hmc_error copy_eoprec_spinor_device(cl_mem in, cl_mem out, usetimer* timer);
	hmc_error convert_to_kappa_format_device(cl_mem inout, const size_t local_work_size, const size_t global_work_size, usetimer* timer);
	hmc_error convert_from_kappa_format_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer);
	hmc_error convert_to_kappa_format_eoprec_device(cl_mem inout, const size_t local_work_size, const size_t global_work_size, usetimer* timer);
	hmc_error convert_from_kappa_format_eoprec_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer);
	hmc_error copy_float_from_device(cl_mem in, hmc_float * out, usetimer* timer);
	hmc_error copy_complex_from_device(cl_mem in, hmc_complex * out, usetimer* timer);
	hmc_error copy_complex_device(cl_mem in, cl_mem out, usetimer* timer);
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
	hmc_error M_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size,  usetimer* timer, usetimer * dslashtimer, usetimer * Mdiagtimer);

	hmc_error bicgstab_device(usetimer * copytimer, usetimer* singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, usetimer * dslashtimer, usetimer * Mdiagtimer, const size_t local_work_size, const size_t global_work_size, int cgmax);
	hmc_error bicgstab_eoprec_device(usetimer * copytimer, usetimer* singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, usetimer * dslashtimer, usetimer * Mdiagtimer, const size_t local_work_size, const size_t global_work_size, int cgmax);
	hmc_error cg_device(usetimer * copytimer, usetimer* singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, usetimer * dslashtimer, usetimer * Mdiagtimer, const size_t local_work_size, const size_t global_work_size, int cgmax);
	hmc_error testing_spinor(inputparameters* parameters, size_t local_size, size_t global_size);

	hmc_error simple_correlator_device(usetimer * copytimer, usetimer* singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, usetimer * solvertimer, usetimer * dslashtimer, usetimer * Mdiagtimer, const size_t ls, const size_t gs, int cgmax);
	hmc_error create_point_source_device(int i, int spacepos, int timepos, const size_t ls, const size_t gs, usetimer * latimer);
	hmc_error create_point_source_eoprec_device(int i, int spacepos, int timepos, const size_t ls, const size_t gs, usetimer * latimer, usetimer * dslashtimer, usetimer * Mdiagtimer);
	hmc_error solver_device(hmc_spinor_field* out, usetimer * copytimer, usetimer * singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, usetimer * dslashtimer, usetimer * Mdiagtimer, usetimer * solvertimer, const size_t ls, const size_t gs, int cgmax);
	hmc_error Aee_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer, usetimer * singletimer, usetimer * dslashtimer, usetimer * Mdiagtimer, usetimer * latimer);
	hmc_error M_inverse_sitediagonal_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer * timer);
	hmc_error M_sitediagonal_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer * timer);
	hmc_error dslash_eoprec_device(cl_mem in, cl_mem out, int evenodd, const size_t local_work_size, const size_t global_work_size, usetimer * timer);
	hmc_error solver_eoprec_device(hmc_spinor_field* out, usetimer * copytimer, usetimer * singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, usetimer * dslashtimer, usetimer * Mdiagtimer, usetimer * solvertimer, const size_t ls, const size_t gs, int cgmax);

	hmc_error perform_benchmark(int cgmax, const size_t ls, const size_t gs, usetimer * copytimer, usetimer * singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, usetimer * solvertimer, usetimer * dslashtimer, usetimer * Mdiagtimer);
	hmc_error finalize_fermions();

#endif
#ifdef _TESTING_
	hmc_error testing(hmc_gaugefield * gaugefield);
#endif
	std::vector<std::string> cl_kernels_file;
	hmc_error finalize();
private:
	int isinit;
	cl_context context;
	cl_command_queue queue;
	cl_program clprogram;
	cl_kernel heatbath_odd;
	cl_kernel heatbath_even;
	cl_kernel overrelax_odd;
	cl_kernel overrelax_even;
	cl_kernel plaquette;
	cl_kernel plaquette_reduction;
	cl_kernel polyakov;
	cl_kernel polyakov_reduction;

	//heatbath variables
	cl_mem clmem_gaugefield;
	cl_mem clmem_rndarray;
	cl_mem clmem_plaq;
	cl_mem clmem_plaq_buf_glob;
	cl_mem clmem_splaq_buf_glob;
	cl_mem clmem_tplaq_buf_glob;
	cl_mem clmem_splaq;
	cl_mem clmem_tplaq;
	cl_mem clmem_polyakov;
	cl_mem clmem_polyakov_buf_glob;
	//!!CP: this is not needed at the moment and since is not copied to the device anywhere!!
	cl_mem clmem_theta_gaugefield;

	//spinorfield and solver variables
#ifdef _FERMIONS_
	cl_mem clmem_kappa;
	cl_mem clmem_theta_fermion;
	cl_mem clmem_mu;
	cl_mem clmem_chem_pot_re;
	cl_mem clmem_chem_pot_im;

	cl_kernel M_diag;
	cl_kernel dslash;
	cl_kernel saxpy;
	cl_kernel saxsbypz;
	cl_kernel scalar_product;
	cl_kernel scalar_product_reduction;
	cl_kernel set_zero_spinorfield;
	cl_kernel global_squarenorm;
	cl_kernel global_squarenorm_reduction;
	cl_kernel ratio;
	cl_kernel product;
	cl_kernel convert_to_kappa_format;
	cl_kernel convert_from_kappa_format;
	cl_kernel convert_to_kappa_format_eoprec;
	cl_kernel convert_from_kappa_format_eoprec;
	cl_kernel create_point_source;
	cl_kernel M_sitediagonal;
	cl_kernel M_inverse_sitediagonal;
	cl_kernel dslash_eoprec;
	cl_kernel saxpy_eoprec;
	cl_kernel saxsbypz_eoprec;
	cl_kernel scalar_product_eoprec;
	cl_kernel set_zero_spinorfield_eoprec;
	cl_kernel global_squarenorm_eoprec;
	cl_kernel create_point_source_eoprec;

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
#endif
	//testing variables
#ifdef _TESTING_
	cl_mem clmem_random_field_int;
	cl_mem clmem_random_field_float;
	cl_mem clmem_random_field_su2;
	cl_mem clmem_heatbath_test_link_in;
	cl_mem clmem_heatbath_test_staple_in;
	cl_mem clmem_heatbath_test_link_out;
	cl_mem clmem_heatbath_test_rnd_array;
	cl_mem clmem_heatbath_test_cter;
	cl_mem clmem_solver_test_spinor_in;
	cl_mem clmem_solver_test_spinor_out;
	cl_mem clmem_solver_test_correlator;
#endif

	/** The number of cores (not PEs) of the device */
	cl_uint max_compute_units;

	/**
	 * Enqueue the given kernel on the device. Local work size will be determined
	 * automatically from device and kernel properties.
	 *
	 * @param kernel The kernel to execute.
	 * @param global_work_size The number of threads to run.
	 *
	 * @todo local work size decision might need ot become less automatic
	 * @todo global work size will also depend on device ...
	 */
	void enqueueKernel(const cl_kernel kernel, const size_t global_work_size);
};

#endif /* _MYOPENCLH_ */
