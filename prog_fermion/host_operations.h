#ifndef _OPERATIONSH_
#define _OPERATIONSH_
#include <iostream>
#include "globaldefs.h"
#include "types.h"
#include "hmcerrs.h"
#include "host_geometry.h"

//local operations
hmc_complex complexconj(hmc_complex *in); 
hmc_complex complexmult(hmc_complex *a, hmc_complex *b); 
hmc_complex complexadd(hmc_complex *a, hmc_complex *b);
hmc_complex complexsubtract(hmc_complex *a, hmc_complex *b); 
hmc_error complexaccumulate(hmc_complex *inout, hmc_complex *incr); 

hmc_complex complexdivide(hmc_complex* numerator, hmc_complex* denominator);

hmc_error adjoin_su3matrix(hmc_su3matrix * mat);
hmc_complex trace_su3matrix(hmc_su3matrix * mat);
hmc_complex det_su3matrix(hmc_su3matrix * U);
hmc_error copy_su3matrix(hmc_su3matrix *out, hmc_su3matrix *in); 
hmc_error copy_staplematrix(hmc_staplematrix *out, hmc_staplematrix *in);
hmc_error unit_su3matrix(hmc_su3matrix * u); 
hmc_error random_su3matrix(hmc_su3matrix * u);
hmc_error zero_su3matrix(hmc_su3matrix * u); 
hmc_error zero_staplematrix(hmc_staplematrix * u);
hmc_error accumulate_su3matrix_prod(hmc_su3matrix *acc, hmc_su3matrix *multiplicator);
hmc_error multiply_su3matrices(hmc_su3matrix *out, hmc_su3matrix *p, hmc_su3matrix *q);
hmc_error multiply_staplematrix(hmc_staplematrix *out, hmc_su3matrix *p, hmc_staplematrix *q);
#ifdef _RECONSTRUCT_TWELVE_
hmc_complex reconstruct_su3(hmc_su3matrix *in, int ncomp);
#endif

//spinor operations, global and local
int spinor_element(int alpha, int color);
int spinor_field_element(int alpha, int color, int nspace, int t);
int spinor_color(int spinor_element);
int spinor_spin(int spinor_element,int color);

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
int eoprec_spinor_field_element(int alpha, int color, int nspace, int t);
int eoprec_spinor_field_element(int alpha, int color, int n_eoprec);
hmc_error convert_from_eoprec(hmc_eoprec_spinor_field* even, hmc_eoprec_spinor_field* odd, hmc_spinor_field* out);
hmc_error convert_to_eoprec(hmc_eoprec_spinor_field* even, hmc_eoprec_spinor_field* odd, hmc_spinor_field* in);
hmc_error convert_to_kappa_format(hmc_eoprec_spinor_field* inout,hmc_float kappa);
hmc_error convert_from_kappa_format(hmc_eoprec_spinor_field* inout,hmc_float kappa);
hmc_float global_squarenorm_eoprec(hmc_eoprec_spinor_field* in);
hmc_complex scalar_product_eoprec(hmc_eoprec_spinor_field* a, hmc_eoprec_spinor_field* b);
hmc_complex scalar_product(hmc_spinor_field* a, hmc_spinor_field* b);
hmc_error get_spinor_from_eoprec_field(hmc_eoprec_spinor_field* in, hmc_spinor* out, int n_eoprec);
hmc_error put_spinor_to_eoprec_field(hmc_spinor* in, hmc_eoprec_spinor_field* out, int n_eoprec);
hmc_error get_spinor_from_field(hmc_spinor_field* in, hmc_spinor* out, int n, int t);
hmc_error put_spinor_to_field(hmc_spinor* in, hmc_spinor_field* out, int n, int t);


//gaugefield and su3 operations, global and local
hmc_error set_gaugefield_cold(hmc_gaugefield * field);
hmc_error set_gaugefield_hot(hmc_gaugefield * field);
hmc_error copy_gaugefield_from_ildg_format(hmc_gaugefield * gaugefield, hmc_float * gaugefield_tmp, int check);
hmc_error copy_gaugefield_to_ildg_format(ildg_gaugefield * dest, hmc_gaugefield * source);

hmc_error copy_to_ocl_format(hmc_ocl_gaugefield* host_gaugefield,hmc_gaugefield* gaugefield);
hmc_error copy_from_ocl_format(hmc_gaugefield* gaugefield,hmc_ocl_gaugefield* host_gaugefield);
int ocl_gaugefield_element(int c, int a, int b, int mu, int spacepos, int t);

hmc_error get_su3matrix(hmc_su3matrix* out, hmc_gaugefield * in, int spacepos, int timepos, int mu); //cl
hmc_error put_su3matrix(hmc_gaugefield * field, hmc_su3matrix * in, int spacepos, int timepos, int mu); //cl

hmc_error adjoin_su3(hmc_gaugefield * in, hmc_gaugefield * out);
hmc_complex global_trace_su3(hmc_gaugefield * field, int mu);

hmc_error accumulate_su3matrices_add(hmc_staplematrix *p, hmc_su3matrix *q);

void reduction (hmc_complex dest[su2_entries], hmc_staplematrix src, const int rand);

void extend (hmc_su3matrix * dest, const int random, hmc_complex src[su2_entries]);

hmc_error project_su3(hmc_su3matrix *U);

#ifdef _USEDOUBLEPREC_
hmc_float const projectioneps = 10.e-12;
#else
hmc_float const projectioneps = 10.e-6;
#endif



#endif
