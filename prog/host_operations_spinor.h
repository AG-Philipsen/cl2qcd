#ifndef _OPERATIONS_SPINORH_
#define _OPERATIONS_SPINORH_
#include <iostream>
#include "globaldefs.h"
#include "types.h"
#include "hmcerrs.h"
#include "host_geometry.h"
#include "host_operations_complex.h"
#include "host_operations_gaugefield.h"
#include <cmath>

hmc_error su3matrix_times_colorvector(hmc_su3matrix* u, hmc_color_vector* in, hmc_color_vector* out);//CP: checked
hmc_error set_local_zero_spinor(hmc_spinor* inout);//CP: checked
hmc_float spinor_squarenorm(hmc_spinor* in);//CP: checked
hmc_error real_multiply_spinor(hmc_spinor* inout, hmc_float factor);//CP: checked, corrected
hmc_error spinors_accumulate(hmc_spinor* inout, hmc_spinor* incr);//CP: checked
hmc_error su3matrix_times_spinor(hmc_su3matrix* u, hmc_spinor* in, hmc_spinor* out);//CP: checked
hmc_error multiply_spinor_gamma0(hmc_spinor* in,hmc_spinor* out);//CP: checked explicitly
hmc_error multiply_spinor_gamma1(hmc_spinor* in,hmc_spinor* out);//CP: checked explicitly
hmc_error multiply_spinor_gamma2(hmc_spinor* in,hmc_spinor* out);//CP: checked explicitly
hmc_error multiply_spinor_gamma3(hmc_spinor* in,hmc_spinor* out);//CP: checked explicitly
hmc_error multiply_spinor_i_factor_gamma5(hmc_spinor* in,hmc_spinor* out, hmc_float factor);//CP: checked explicitly
hmc_error spinprojectproduct_gamma0(hmc_su3matrix* u, hmc_spinor* spin, hmc_float sign);//CP: checked explicitly
hmc_error spinprojectproduct_gamma1(hmc_su3matrix* u, hmc_spinor* spin, hmc_float sign);//CP: checked explicitly
hmc_error spinprojectproduct_gamma2(hmc_su3matrix* u, hmc_spinor* spin, hmc_float sign);//CP: checked explicitly
hmc_error spinprojectproduct_gamma3(hmc_su3matrix* u, hmc_spinor* spin, hmc_float sign);//CP: checked explicitly
hmc_float spinor_squarenorm(hmc_spinor* in);//CP: checked explicitly
hmc_error spinors_accumulate(hmc_spinor* inout, hmc_spinor* incr);//CP: checked explicitly
void spinor_apply_bc(hmc_spinor * in, hmc_float theta);//CP: checked explicitly

#endif
