#ifndef _OPERATIONS_SPINORH_
#define _OPERATIONS_SPINORH_
#include <iostream>
#include "globaldefs.h"
#include "types.h"
#include "hmcerrs.h"
#include "host_geometry.h"
#include "host_operations_complex.h"
#include "host_operations_gaugefield.h"
#include "host_solver.h"
#include <cmath>

//spinor operations, global and local
hmc_error set_zero_spinor(hmc_spinor_field * field);
hmc_error set_zero_eoprec_spinor(hmc_eoprec_spinor_field * field);
hmc_float local_squarenorm(hmc_spinor_field *field, int spacepos, int timepos);
hmc_float global_squarenorm(hmc_spinor_field *field);
hmc_error fill_with_one(hmc_spinor_field *field, int spacepos, int timepos, int alpha, int color);

hmc_error su3matrix_times_colorvector(hmc_su3matrix* u, hmc_color_vector* in, hmc_color_vector* out);
hmc_error set_local_zero_spinor(hmc_spinor* inout);
hmc_float spinor_squarenorm(hmc_spinor* in);
hmc_error real_multiply_spinor(hmc_spinor* inout, hmc_float factor);

hmc_error spinprojectproduct_gamma0(hmc_su3matrix* u, hmc_spinor* spin, hmc_float sign);
hmc_error spinprojectproduct_gamma1(hmc_su3matrix* u, hmc_spinor* spin, hmc_float sign);
hmc_error spinprojectproduct_gamma2(hmc_su3matrix* u, hmc_spinor* spin, hmc_float sign);
hmc_error spinprojectproduct_gamma3(hmc_su3matrix* u, hmc_spinor* spin, hmc_float sign);

hmc_error spinors_accumulate(hmc_spinor* inout, hmc_spinor* incr);
hmc_error multiply_spinor_factor_gamma5(hmc_spinor* in,hmc_spinor* out, hmc_float factor);
hmc_error multiply_spinor_gamma0(hmc_spinor* in,hmc_spinor* out);
hmc_error multiply_spinor_gamma1(hmc_spinor* in,hmc_spinor* out);
hmc_error multiply_spinor_gamma2(hmc_spinor* in,hmc_spinor* out);
hmc_error multiply_spinor_gamma3(hmc_spinor* in,hmc_spinor* out);
hmc_error su3matrix_times_spinor(hmc_su3matrix* u, hmc_spinor* in, hmc_spinor* out);


//eoprec operations
hmc_error convert_from_eoprec(hmc_eoprec_spinor_field* even, hmc_eoprec_spinor_field* odd, hmc_spinor_field* out);
hmc_error convert_to_eoprec(hmc_eoprec_spinor_field* even, hmc_eoprec_spinor_field* odd, hmc_spinor_field* in);
hmc_error convert_to_kappa_format_eoprec(hmc_eoprec_spinor_field* inout,hmc_float kappa);
hmc_error convert_from_kappa_format_eoprec(hmc_eoprec_spinor_field* in, hmc_eoprec_spinor_field * out,hmc_float kappa);
hmc_error convert_to_kappa_format(hmc_spinor_field* inout,hmc_float kappa);
hmc_error convert_from_kappa_format(hmc_spinor_field* in, hmc_spinor_field * out,hmc_float kappa);
hmc_float global_squarenorm_eoprec(hmc_eoprec_spinor_field* in);
hmc_complex scalar_product_eoprec(hmc_eoprec_spinor_field* a, hmc_eoprec_spinor_field* b);
hmc_complex scalar_product(hmc_spinor_field* a, hmc_spinor_field* b);
hmc_error get_spinor_from_eoprec_field(hmc_eoprec_spinor_field* in, hmc_spinor* out, int n_eoprec);
hmc_error put_spinor_to_eoprec_field(hmc_spinor* in, hmc_eoprec_spinor_field* out, int n_eoprec);
hmc_error get_spinor_from_field(hmc_spinor_field* in, hmc_spinor* out, int n, int t);
hmc_error put_spinor_to_field(hmc_spinor* in, hmc_spinor_field* out, int n, int t);

//new stuff

void copy_spinor(hmc_complex * in, hmc_complex * out);
void copy_spinor_eoprec(hmc_complex * in, hmc_complex * out);

// -alpha*x + y
//CP: defined with a minus!!!
void saxpy(hmc_spinor_field * x, hmc_spinor_field * y, hmc_complex * alpha, hmc_spinor_field * out);
void saxpy_eoprec(hmc_eoprec_spinor_field * x, hmc_eoprec_spinor_field * y, hmc_complex * alpha, hmc_eoprec_spinor_field * out);

//alpha*x + beta*y + z
void saxsbypz(hmc_spinor_field * x, hmc_spinor_field * y,  hmc_spinor_field * z, hmc_complex * alpha, hmc_complex * beta, hmc_spinor_field * out);
void saxsbypz_eoprec(hmc_eoprec_spinor_field * x, hmc_eoprec_spinor_field * y,  hmc_eoprec_spinor_field * z, hmc_complex * alpha, hmc_complex * beta, hmc_eoprec_spinor_field * out);

hmc_error create_point_source(hmc_spinor_field* b, int i, int spacepos, int timepos, hmc_float kappa, hmc_float mu, hmc_gaugefield* gaugefield);
hmc_error create_point_source_eoprec(hmc_eoprec_spinor_field* be,hmc_eoprec_spinor_field* bo,int i,int spacepos,int timepos,hmc_float kappa, hmc_float mu, hmc_float theta, hmc_gaugefield* gaugefield);

void spinor_apply_bc(hmc_spinor * in, hmc_float theta);

#endif
