#ifndef _MYOPENCLH_
#define _MYOPENCLH_

#include <cstdlib>
#include <vector>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <CL/cl.h>

#include "host_operations.h"
#include "globaldefs.h"
#include "hmcerrs.h"
#include "types.h"
#include "host_use_timer.h"
#include "host_testing.h"
#include "host_random.h"

//give a list of all kernel-files
std::vector<std::string> const cl_kernels_file = {"opencl_header.cl", "opencl_geometry.cl","opencl_operations.cl",   "opencl_random.cl", "opencl_update_heatbath.cl", "opencl_solver.cl", "opencl_gaugeobservables.cl",  "opencl_fermionobservables.cl", "opencl_testing.cl"};

class opencl {
 public:
  opencl(cl_device_type wanted, usetimer* timer){init(wanted, timer);};
  ~opencl(){finalize();};
  hmc_error init(cl_device_type wanted_device_type, usetimer* timer);
  hmc_error copy_gaugefield_to_device(hmc_gaugefield* host_gaugefield,  usetimer* timer);
  hmc_error copy_rndarray_to_device(hmc_rndarray host_rndarray,  usetimer* timer);
  hmc_error copy_rndarray_from_device(hmc_rndarray rndarray, usetimer* timer);
  hmc_error get_gaugefield_from_device(hmc_gaugefield* host_gaugefield,  usetimer* timer);
  hmc_error run_heatbath(double beta, const size_t local_work_size, const size_t global_work_size,  usetimer* timer);
  hmc_error run_overrelax(double beta, const size_t local_work_size, const size_t global_work_size,  usetimer* timer);
  hmc_error gaugeobservables(const size_t local_work_size, const size_t global_work_size, hmc_float * plaq, hmc_float * tplaq, hmc_float * splaq, hmc_complex * pol, usetimer* timer1, usetimer* timer2);
  hmc_error testing(hmc_gaugefield * gaugefield);
  hmc_error finalize();
 private:
  int isinit;
  cl_context context;
  cl_command_queue queue;
  cl_program clprogram;
  cl_mem clmem_gaugefield;
  cl_mem clmem_rndarray;
  //CP: unused at the moment
  cl_mem clmem_spinorfield;
  cl_kernel heatbath_odd;
  cl_kernel heatbath_even;
  cl_kernel overrelax_odd;
  cl_kernel overrelax_even;
  cl_kernel plaquette;
  cl_kernel polyakov;
  cl_mem clmem_random_field_int;
  cl_mem clmem_random_field_float;
  cl_mem clmem_random_field_su2;
  cl_mem clmem_plaq;
  cl_mem clmem_splaq;
  cl_mem clmem_tplaq;
  cl_mem clmem_polyakov;
  cl_mem clmem_heatbath_test_link_in;
  cl_mem clmem_heatbath_test_staple_in;
  cl_mem clmem_heatbath_test_link_out;
  cl_mem clmem_heatbath_test_rnd_array;
  cl_mem clmem_heatbath_test_cter;
  cl_mem clmem_solver_test_spinor_in;
  cl_mem clmem_solver_test_spinor_out;
  cl_mem clmem_solver_test_correlator;
};

#endif
