#ifndef _MYOPENCLH_
#define _MYOPENCLH_

#include <cstdlib>
#include <vector>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <CL/cl.h>

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

class opencl {
 public:
  opencl(cl_device_type wanted, const size_t ls, const size_t gs, usetimer* timer){init(wanted, ls, gs, timer);};
  ~opencl(){finalize();};
  hmc_error init(cl_device_type wanted_device_type, const size_t local_work_size, const size_t global_work_size, usetimer* timer);
  hmc_error copy_gaugefield_to_device(hmc_gaugefield* host_gaugefield,  usetimer* timer);
  hmc_error copy_rndarray_to_device(hmc_rndarray host_rndarray,  usetimer* timer);
  hmc_error copy_rndarray_from_device(hmc_rndarray rndarray, usetimer* timer);
  hmc_error get_gaugefield_from_device(hmc_gaugefield* host_gaugefield,  usetimer* timer);
  hmc_error run_heatbath(hmc_float beta, const size_t local_work_size, const size_t global_work_size,  usetimer* timer);
  hmc_error run_overrelax(hmc_float beta, const size_t local_work_size, const size_t global_work_size,  usetimer* timer);
  hmc_error gaugeobservables(const size_t local_work_size, const size_t global_work_size, hmc_float * plaq, hmc_float * tplaq, hmc_float * splaq, hmc_complex * pol, usetimer* timer1, usetimer* timer2);
#ifdef _FERMIONS_
	hmc_error init_solver_variables(inputparameters* parameters, const size_t local_work_size, const size_t global_work_size, usetimer* timer);
	hmc_error copy_spinorfield_to_device(hmc_spinor_field* host_spinorfield, usetimer* timer);
	hmc_error copy_source_to_device(hmc_spinor_field* host_spinorfield, usetimer* timer);
	hmc_error get_spinorfield_from_device(hmc_spinor_field* host_spinorfield,  usetimer* timer);
	hmc_error copy_spinor_device(cl_mem in, cl_mem out, usetimer* timer);
	hmc_error convert_to_kappa_format_device(cl_mem inout, const size_t local_work_size, const size_t global_work_size, usetimer* timer);
	hmc_error convert_from_kappa_format_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer);
	hmc_error convert_to_kappa_format_eoprec_device(cl_mem inout, const size_t local_work_size, const size_t global_work_size, usetimer* timer);
	hmc_error convert_from_kappa_format_eoprec_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer);
	hmc_error copy_float_from_device(cl_mem in, hmc_float * out, usetimer* timer);
	hmc_error copy_complex_from_device(cl_mem in, hmc_complex * out, usetimer* timer);
	hmc_error copy_complex_device(cl_mem in, cl_mem out, usetimer* timer);
	hmc_error set_complex_to_scalar_product_device(cl_mem a, cl_mem b, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer);
	hmc_error set_complex_to_ratio_device(cl_mem a, cl_mem b, cl_mem out, usetimer* timer);
	hmc_error set_complex_to_product_device(cl_mem a, cl_mem b, cl_mem out, usetimer* timer);
	hmc_error set_float_to_global_squarenorm_device(cl_mem a, cl_mem out, const size_t local_work_size, const size_t global_work_size, usetimer* timer);
	hmc_error set_zero_spinorfield_device(cl_mem x, const size_t local_work_size, const size_t global_work_size, usetimer* timer);
	hmc_error saxpy_device(cl_mem x, cl_mem y, cl_mem alpha, cl_mem out, const size_t local_work_size, const size_t global_work_size,  usetimer* timer);
	hmc_error saxsbypz_device(cl_mem x, cl_mem y, cl_mem z, cl_mem alpha, cl_mem beta, cl_mem out, const size_t local_work_size, const size_t global_work_size,  usetimer* timer);
	hmc_error M_device(cl_mem in, cl_mem out, const size_t local_work_size, const size_t global_work_size,  usetimer* timer);

	hmc_error bicgstab_device(usetimer * copytimer, usetimer* singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, const size_t local_work_size, const size_t global_work_size, int cgmax);
	hmc_error cg_device(usetimer * copytimer, usetimer* singletimer, usetimer * Mtimer, usetimer * scalarprodtimer, usetimer * latimer, const size_t local_work_size, const size_t global_work_size, int cgmax);
	hmc_error testing_spinor(inputparameters* parameters, size_t local_size, size_t global_size);
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
	//!!CP: buffer for the local scalar_product
	cl_mem clmem_scalar_product_buf_loc;// [local_work_size];
	//!!CP: buffer for the global scalar_product
	cl_mem clmem_scalar_product_buf_glob;//[num_group];
	//!!CP: buffer for the local global_squarenorm_product
	cl_mem clmem_global_squarenorm_buf_loc;// [local_work_size];
	//!!CP: buffer for the global global_squarenorm_product
	cl_mem clmem_global_squarenorm_buf_glob;//[num_group];
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

};

#endif
