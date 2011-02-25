#include "host_operations_spinor.h"

hmc_error su3matrix_times_colorvector(hmc_su3matrix* u, hmc_color_vector* in, hmc_color_vector* out){
#ifdef _RECONSTRUCT_TWELVE_
  for(int a=0; a<NC-1; a++) {
    out[a] = hmc_complex_zero;
    for(int b=0; b<NC; b++) {
      hmc_complex tmp = complexmult(&((*u)[a+(NC-1)*b]),&(in[b]));
      complexaccumulate(&out[a],&tmp);
    }
  }
  out[2] = hmc_complex_zero;
  for(int b=0; b<NC; b++) {
    hmc_complex rec = reconstruct_su3(u,b);
    hmc_complex tmp = complexmult(&rec,&(in[b]));
    complexaccumulate(&out[2],&tmp);
  }
#else
  for(int a=0; a<NC; a++) {
    out[a] = hmc_complex_zero;
    for(int b=0; b<NC; b++) {
      hmc_complex tmp2 = (*u)[a][b];
      hmc_complex tmp = complexmult(&tmp2,&(in[b]));
      complexaccumulate(&(out[a]),&tmp);
    }
  }
#endif
  return HMC_SUCCESS;
}

hmc_error set_zero_spinor(hmc_spinor_field *field) {
  for (int n=0; n<SPINORFIELDSIZE; n++) {
	field[n].re=0;
	field[n].im=0;
  }
  return HMC_SUCCESS;
}

hmc_error set_zero_eoprec_spinor(hmc_eoprec_spinor_field *field) {
  for (int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
	field[n].re=0;
	field[n].im=0;
  }
  return HMC_SUCCESS;
}

hmc_float local_squarenorm(hmc_spinor_field *field, int spacepos, int timepos) {
  hmc_float sum=0;
  for (int a=0; a<NSPIN; a++) {
    for(int c=0; c<NC; c++) {
      hmc_float dummy_re = field[spinor_field_element(a,c,spacepos,timepos)].re;
      hmc_float dummy_im = field[spinor_field_element(a,c,spacepos,timepos)].im;
      sum += dummy_re*dummy_re + dummy_im*dummy_im;
    }
  }
  return sum;
}

hmc_float global_squarenorm(hmc_spinor_field *field) {
  hmc_float sum=0;
  for (int t=0; t<NTIME; t++) {
    for (int n=0; n<VOLSPACE; n++) {
      sum += local_squarenorm(field,n,t);
    }
  }
  return sum;
}

hmc_error fill_with_one(hmc_spinor_field *field, int spacepos, int timepos, int alpha, int color){
  field[spinor_field_element(alpha,color,spacepos,timepos)].re = hmc_one_f;
  field[spinor_field_element(alpha,color,spacepos,timepos)].im = 0;
  return HMC_SUCCESS;
}

hmc_error convert_to_eoprec(hmc_eoprec_spinor_field* even, hmc_eoprec_spinor_field* odd, hmc_spinor_field* in){
  int spacepos;
  int timepos;
  for(int n=0; n<VOL4D/2; n++) {
    for(int alpha=0; alpha<NSPIN; alpha++) {
      for(int color=0; color<NC; color++) {
	get_even_site(n, &spacepos, &timepos);
	even[eoprec_spinor_field_element(alpha,color,n)] = in[spinor_field_element(alpha,color,spacepos,timepos)];
	get_odd_site(n, &spacepos, &timepos);
	odd[eoprec_spinor_field_element(alpha,color,n)] = in[spinor_field_element(alpha,color,spacepos,timepos)];
      }
    }
  }
  return HMC_SUCCESS;
}

hmc_error convert_from_eoprec(hmc_eoprec_spinor_field* even, hmc_eoprec_spinor_field* odd, hmc_spinor_field* out){
  int spacepos, timepos;
  for(int n=0; n<VOL4D/2; n++) {
    for(int alpha=0; alpha<NSPIN; alpha++) {
      for(int color=0; color<NC; color++) {
        get_even_site(n, &spacepos, &timepos);
        out[spinor_field_element(alpha,color,spacepos,timepos)] = even[eoprec_spinor_field_element(alpha,color,n)];
	get_odd_site(n, &spacepos, &timepos);
	out[spinor_field_element(alpha,color,spacepos,timepos)] = odd[eoprec_spinor_field_element(alpha,color,n)];
      }
    }
  }
  return HMC_SUCCESS;
}

hmc_error convert_to_kappa_format(hmc_spinor_field* inout,hmc_float kappa){
  for(int n=0; n<SPINORFIELDSIZE; n++) {
    inout[n].re *= sqrt(2.*kappa);
    inout[n].im *= sqrt(2.*kappa);
  }
  return HMC_SUCCESS;
}

hmc_error convert_to_kappa_format_eoprec(hmc_eoprec_spinor_field* inout,hmc_float kappa){
  for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
    inout[n].re *= sqrt(2.*kappa);
    inout[n].im *= sqrt(2.*kappa);
  }
  return HMC_SUCCESS;
}

hmc_error convert_from_kappa_format(hmc_spinor_field* in, hmc_spinor_field * out,hmc_float kappa){
  for(int n=0; n<SPINORFIELDSIZE; n++) {
    out[n].re = (in[n].re)/sqrt(2.*kappa);
    out[n].im = (in[n].im)/sqrt(2.*kappa);
  }
  return HMC_SUCCESS;
}

hmc_error convert_from_kappa_format_eoprec(hmc_eoprec_spinor_field* in, hmc_eoprec_spinor_field * out, hmc_float kappa){
  for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
    out[n].re = (in[n].re)/sqrt(2.*kappa);
    out[n].im = (in[n].im)/sqrt(2.*kappa);
  }
  return HMC_SUCCESS;
}

hmc_float global_squarenorm_eoprec(hmc_eoprec_spinor_field *in) {
  hmc_float sum=0;
  for (int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
    sum += in[n].re*in[n].re + in[n].im*in[n].im;
  }
  return sum;
}

hmc_complex scalar_product_eoprec(hmc_eoprec_spinor_field* a, hmc_eoprec_spinor_field* b){
  // (a,b) = sum_k conj(a_k)*b_k
  hmc_complex res;
  res.re=0;
  res.im=0;
  for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
    res.re += a[n].re*b[n].re + a[n].im*b[n].im;
    res.im += a[n].re*b[n].im - a[n].im*b[n].re;
  }
  return res;
}

hmc_complex scalar_product(hmc_spinor_field* a, hmc_spinor_field* b){
  // (a,b) = sum_k conj(a_k)*b_k
  hmc_complex res;
  res.re=0;
  res.im=0;
  for(int n=0; n<SPINORFIELDSIZE; n++) {
    res.re += a[n].re*b[n].re + a[n].im*b[n].im;
    res.im += a[n].re*b[n].im - a[n].im*b[n].re;
  }
  return res;
}

hmc_error get_spinor_from_eoprec_field(hmc_eoprec_spinor_field* in, hmc_spinor* out, int n_eoprec){
  for(int alpha=0; alpha<NSPIN; alpha++) {
    for(int color=0; color<NC; color++) {
      out[spinor_element(alpha,color)] = in[eoprec_spinor_field_element(alpha,color,n_eoprec)];
    }
  }
  return HMC_SUCCESS;
}

hmc_error put_spinor_to_eoprec_field(hmc_spinor* in, hmc_eoprec_spinor_field* out, int n_eoprec){
  for(int alpha=0; alpha<NSPIN; alpha++) {
    for(int color=0; color<NC; color++) {
      out[eoprec_spinor_field_element(alpha,color,n_eoprec)]=in[spinor_element(alpha,color)];
    }
  }
  return HMC_SUCCESS;
}

hmc_error get_spinor_from_field(hmc_spinor_field* in, hmc_spinor* out, int n, int t){
  for(int alpha=0; alpha<NSPIN; alpha++) {
    for(int color=0; color<NC; color++) {
      out[spinor_element(alpha,color)] = in[spinor_field_element(alpha,color,n,t)];
    }
  }
  return HMC_SUCCESS;
}

hmc_error put_spinor_to_field(hmc_spinor* in, hmc_spinor_field* out, int n, int t){
  for(int alpha=0; alpha<NSPIN; alpha++) {
    for(int color=0; color<NC; color++) {
      out[spinor_field_element(alpha,color,n,t)]=in[spinor_element(alpha,color)];
    }
  }
  return HMC_SUCCESS;
}

hmc_error set_local_zero_spinor(hmc_spinor* inout){
  for(int j=0; j<SPINORSIZE; j++) {
    inout[j].re = 0;
    inout[j].im = 0;
  }
  return HMC_SUCCESS;
}

hmc_float spinor_squarenorm(hmc_spinor* in){
  hmc_float res=0;
  for(int j=0; j<SPINORSIZE; j++) {
    hmc_complex tmp = complexconj(&in[j]);
    hmc_complex incr= complexmult(&tmp,&in[j]);
    res+=incr.re;
  }
  return res;
}

hmc_error real_multiply_spinor(hmc_spinor* inout, hmc_float factor){
  for(int j=0; j<SPINORSIZE; j++) inout[j].re *=factor;
  return HMC_SUCCESS;
}

hmc_error spinprojectproduct_gamma0(hmc_su3matrix* u, hmc_spinor* spin,hmc_float sign){
  //out = u*(1+sign*gamma0)*in
  hmc_color_vector vec1[NC];
  hmc_color_vector vec2[NC];
  hmc_color_vector uvec1[NC];
  hmc_color_vector uvec2[NC];
  for(int c=0; c<NC; c++) {
    vec1[c].re = spin[spinor_element(0,c)].re - sign*spin[spinor_element(3,c)].im;
    vec1[c].im = spin[spinor_element(0,c)].im + sign*spin[spinor_element(3,c)].re;
    vec2[c].re = spin[spinor_element(1,c)].re - sign*spin[spinor_element(2,c)].im;
    vec2[c].im = spin[spinor_element(1,c)].im + sign*spin[spinor_element(2,c)].re;
  }
  su3matrix_times_colorvector(u,vec1,uvec1);
  su3matrix_times_colorvector(u,vec2,uvec2);
  for(int c=0; c<NC; c++) {
    spin[spinor_element(0,c)].re = uvec1[c].re;
    spin[spinor_element(0,c)].im = uvec1[c].im;
    spin[spinor_element(1,c)].re = uvec2[c].re;
    spin[spinor_element(1,c)].im = uvec2[c].im;
    spin[spinor_element(2,c)].re = sign*uvec2[c].im;
    spin[spinor_element(2,c)].im = -sign*uvec2[c].re;
    spin[spinor_element(3,c)].re = sign*uvec1[c].im;
    spin[spinor_element(3,c)].im = -sign*uvec1[c].re;
  }
  return HMC_SUCCESS;
}

hmc_error spinprojectproduct_gamma1(hmc_su3matrix* u, hmc_spinor* spin, hmc_float sign){
  //out = u*(1+sign*gamma1)*in
  hmc_color_vector vec1[NC];
  hmc_color_vector vec2[NC];
  hmc_color_vector uvec1[NC];
  hmc_color_vector uvec2[NC];
  for(int c=0; c<NC; c++) {
    vec1[c].re = spin[spinor_element(0,c)].re - sign*spin[spinor_element(3,c)].re;
    vec1[c].im = spin[spinor_element(0,c)].im - sign*spin[spinor_element(3,c)].im;
    vec2[c].re = spin[spinor_element(1,c)].re + sign*spin[spinor_element(2,c)].re;
    vec2[c].im = spin[spinor_element(1,c)].im + sign*spin[spinor_element(2,c)].im;
  }
  su3matrix_times_colorvector(u,vec1,uvec1);
  su3matrix_times_colorvector(u,vec2,uvec2);
  for(int c=0; c<NC; c++) {
    spin[spinor_element(0,c)].re = uvec1[c].re;
    spin[spinor_element(0,c)].im = uvec1[c].im;
    spin[spinor_element(1,c)].re = uvec2[c].re;
    spin[spinor_element(1,c)].im = uvec2[c].im;
    spin[spinor_element(2,c)].re = sign*uvec2[c].re;
    spin[spinor_element(2,c)].im = sign*uvec2[c].im;
    spin[spinor_element(3,c)].re = -sign*uvec1[c].re;
    spin[spinor_element(3,c)].im = -sign*uvec1[c].im;
  }
  return HMC_SUCCESS;
}

hmc_error spinprojectproduct_gamma2(hmc_su3matrix* u, hmc_spinor* spin, hmc_float sign){
  //out = u*(1+sign*gamma2)*in
  hmc_color_vector vec1[NC];
  hmc_color_vector vec2[NC];
  hmc_color_vector uvec1[NC];
  hmc_color_vector uvec2[NC];
  for(int c=0; c<NC; c++) {
    vec1[c].re = spin[spinor_element(0,c)].re - sign*spin[spinor_element(2,c)].im;
    vec1[c].im = spin[spinor_element(0,c)].im + sign*spin[spinor_element(2,c)].re;
    vec2[c].re = spin[spinor_element(1,c)].re + sign*spin[spinor_element(3,c)].im;
    vec2[c].im = spin[spinor_element(1,c)].im - sign*spin[spinor_element(3,c)].re;
  }
  su3matrix_times_colorvector(u,vec1,uvec1);
  su3matrix_times_colorvector(u,vec2,uvec2);
  for(int c=0; c<NC; c++) {
    spin[spinor_element(0,c)].re = uvec1[c].re;
    spin[spinor_element(0,c)].im = uvec1[c].im;
    spin[spinor_element(1,c)].re = uvec2[c].re;
    spin[spinor_element(1,c)].im = uvec2[c].im;
    spin[spinor_element(2,c)].re = sign*uvec1[c].im;
    spin[spinor_element(2,c)].im = -sign*uvec1[c].re;
    spin[spinor_element(3,c)].re = -sign*uvec2[c].im;
    spin[spinor_element(3,c)].im = sign*uvec2[c].re;
  }
  return HMC_SUCCESS;
}

hmc_error spinprojectproduct_gamma3(hmc_su3matrix* u, hmc_spinor* spin, hmc_float sign){
  //out = u*(1+sign*gamma3)*in
  hmc_color_vector vec1[NC];
  hmc_color_vector vec2[NC];
  hmc_color_vector uvec1[NC];
  hmc_color_vector uvec2[NC];
  for(int c=0; c<NC; c++) {
    vec1[c].re = spin[spinor_element(0,c)].re + sign*spin[spinor_element(2,c)].re;
    vec1[c].im = spin[spinor_element(0,c)].im + sign*spin[spinor_element(2,c)].im;
    vec2[c].re = spin[spinor_element(1,c)].re + sign*spin[spinor_element(3,c)].re;
    vec2[c].im = spin[spinor_element(1,c)].im + sign*spin[spinor_element(3,c)].im;
  }
  su3matrix_times_colorvector(u,vec1,uvec1);
  su3matrix_times_colorvector(u,vec2,uvec2);
  for(int c=0; c<NC; c++) {
    spin[spinor_element(0,c)].re = uvec1[c].re;
    spin[spinor_element(0,c)].im = uvec1[c].im;
    spin[spinor_element(1,c)].re = uvec2[c].re;
    spin[spinor_element(1,c)].im = uvec2[c].im;
    spin[spinor_element(2,c)].re = sign*uvec1[c].re;
    spin[spinor_element(2,c)].im = sign*uvec1[c].im;
    spin[spinor_element(3,c)].re = sign*uvec2[c].re;
    spin[spinor_element(3,c)].im = sign*uvec2[c].im;
  }
  return HMC_SUCCESS;
}

hmc_error spinors_accumulate(hmc_spinor* inout, hmc_spinor* incr){
  for(int j=0; j<SPINORSIZE; j++) {
    inout[j].re += incr[j].re;
    inout[j].im += incr[j].im;
  }
  return HMC_SUCCESS;
}

hmc_error multiply_spinor_factor_gamma5(hmc_spinor* in, hmc_spinor* out, hmc_float factor){
  for(int c=0; c<NC; c++) {
    out[spinor_element(0,c)].re = -factor*in[spinor_element(0,c)].im;
    out[spinor_element(1,c)].re = -factor*in[spinor_element(1,c)].im;
    out[spinor_element(2,c)].re = factor*in[spinor_element(2,c)].im;
    out[spinor_element(3,c)].re = factor*in[spinor_element(3,c)].im;
    out[spinor_element(0,c)].im = factor*in[spinor_element(0,c)].re;
    out[spinor_element(1,c)].im = factor*in[spinor_element(1,c)].re;
    out[spinor_element(2,c)].im = -factor*in[spinor_element(2,c)].re;
    out[spinor_element(3,c)].im = -factor*in[spinor_element(3,c)].re;
  }
  return HMC_SUCCESS;
}

hmc_error multiply_spinor_gamma0(hmc_spinor* in,hmc_spinor* out){
  for(int c=0; c<NC; c++) {
    out[spinor_element(0,c)].re = -in[spinor_element(3,c)].im;
    out[spinor_element(1,c)].re = -in[spinor_element(2,c)].im;
    out[spinor_element(2,c)].re = in[spinor_element(1,c)].im;
    out[spinor_element(3,c)].re = in[spinor_element(0,c)].im;
    out[spinor_element(0,c)].im = in[spinor_element(3,c)].re;
    out[spinor_element(1,c)].im = in[spinor_element(2,c)].re;
    out[spinor_element(2,c)].im = -in[spinor_element(1,c)].re;
    out[spinor_element(3,c)].im = -in[spinor_element(0,c)].re;
  }
  return HMC_SUCCESS;
}
hmc_error multiply_spinor_gamma1(hmc_spinor* in,hmc_spinor* out){
  for(int c=0; c<NC; c++) {
    out[spinor_element(0,c)].re = -in[spinor_element(3,c)].re;
    out[spinor_element(1,c)].re = in[spinor_element(2,c)].re;
    out[spinor_element(2,c)].re = in[spinor_element(1,c)].re;
    out[spinor_element(3,c)].re = -in[spinor_element(0,c)].re;
    out[spinor_element(0,c)].im = -in[spinor_element(3,c)].im;
    out[spinor_element(1,c)].im = in[spinor_element(2,c)].im;
    out[spinor_element(2,c)].im = in[spinor_element(1,c)].im;
    out[spinor_element(3,c)].im = -in[spinor_element(0,c)].im;
  }
  return HMC_SUCCESS;
}
hmc_error multiply_spinor_gamma2(hmc_spinor* in,hmc_spinor* out){
  for(int c=0; c<NC; c++) {
    out[spinor_element(0,c)].re = -in[spinor_element(2,c)].im;
    out[spinor_element(1,c)].re = in[spinor_element(3,c)].im;
    out[spinor_element(2,c)].re = in[spinor_element(0,c)].im;
    out[spinor_element(3,c)].re = -in[spinor_element(1,c)].im;
    out[spinor_element(0,c)].im = in[spinor_element(2,c)].re;
    out[spinor_element(1,c)].im = -in[spinor_element(3,c)].re;
    out[spinor_element(2,c)].im = -in[spinor_element(0,c)].re;
    out[spinor_element(3,c)].im = in[spinor_element(1,c)].re;
  }
  return HMC_SUCCESS;
}
hmc_error multiply_spinor_gamma3(hmc_spinor* in,hmc_spinor* out){
  for(int c=0; c<NC; c++) {
    out[spinor_element(0,c)].re = in[spinor_element(2,c)].re;
    out[spinor_element(1,c)].re = in[spinor_element(3,c)].re;
    out[spinor_element(2,c)].re = in[spinor_element(0,c)].re;
    out[spinor_element(3,c)].re = in[spinor_element(1,c)].re;
    out[spinor_element(0,c)].im = in[spinor_element(2,c)].im;
    out[spinor_element(1,c)].im = in[spinor_element(3,c)].im;
    out[spinor_element(2,c)].im = in[spinor_element(0,c)].im;
    out[spinor_element(3,c)].im = in[spinor_element(1,c)].im;
  }
  return HMC_SUCCESS;
}

hmc_error su3matrix_times_spinor(hmc_su3matrix* u, hmc_spinor* in, hmc_spinor* out){
  for (int alpha=0; alpha<NSPIN; alpha++) {
    hmc_color_vector vec_in[NC];
    hmc_color_vector vec_out[NC];
    for(int c=0; c<NC; c++) {
      vec_in[c].re = in[spinor_element(alpha,c)].re;
      vec_in[c].im = in[spinor_element(alpha,c)].im;
    }
    su3matrix_times_colorvector(u, vec_in, vec_out);
    for(int c=0; c<NC; c++) {
      in[spinor_element(alpha,c)].re = vec_out[c].re;
      in[spinor_element(alpha,c)].im = vec_out[c].im;
    }
  }
  return HMC_SUCCESS;
}

void copy_spinor(hmc_complex * in, hmc_complex * out){
  for (int n=0; n<SPINORFIELDSIZE; n++) {
    out[n].re = in[n].re;
    out[n].im = in[n].im;
  }
  return;
}

void copy_spinor_eoprec(hmc_complex * in, hmc_complex * out){
  for (int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
    out[n].re = in[n].re;
    out[n].im = in[n].im;
  }
  return;
}

// -alpha*x + y
//CP: defined with a minus!!!
void saxpy(hmc_spinor_field * x, hmc_spinor_field * y, hmc_complex * alpha, hmc_spinor_field * out){
  for (int n=0; n<SPINORFIELDSIZE; n++) {
    hmc_complex tmp1 = complexmult(alpha,&x[n]);
    ((out)[n]).re = -(tmp1).re + y[n].re;
    ((out)[n]).im = -(tmp1).im + y[n].im;
  }
  return;
}

void saxpy_eoprec(hmc_eoprec_spinor_field * x, hmc_eoprec_spinor_field * y, hmc_complex * alpha, hmc_eoprec_spinor_field * out){
  for (int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
    hmc_complex tmp1 = complexmult(alpha,&x[n]);
    ((out)[n]).re = -(tmp1).re + y[n].re;
    ((out)[n]).im = -(tmp1).im + y[n].im;
  }
  return;
}

//alpha*x + beta*y + z
void saxsbypz(hmc_spinor_field * x, hmc_spinor_field * y,  hmc_spinor_field * z, hmc_complex * alpha, hmc_complex * beta, hmc_spinor_field * out){
  for (int n=0; n<SPINORFIELDSIZE; n++) {
    hmc_complex tmp1 = complexmult(alpha,&x[n]);
    hmc_complex tmp2 = complexmult(beta,&y[n]);
    ((out)[n]).re = (tmp1).re + (tmp2).re + z[n].re;
    ((out)[n]).im = (tmp1).im + (tmp2).im + z[n].im;
  }
  return;
}

void saxsbypz_eoprec(hmc_eoprec_spinor_field * x, hmc_eoprec_spinor_field * y,  hmc_eoprec_spinor_field * z, hmc_complex * alpha, hmc_complex * beta, hmc_eoprec_spinor_field * out){
  for (int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
    hmc_complex tmp1 = complexmult(alpha,&x[n]);
    hmc_complex tmp2 = complexmult(beta,&y[n]);
    ((out)[n]).re = (tmp1).re + (tmp2).re + z[n].re;
    ((out)[n]).im = (tmp1).im + (tmp2).im + z[n].im;
  }
  return;
}

hmc_error create_point_source(hmc_spinor_field* b, int i, int spacepos, int timepos, hmc_float kappa, hmc_float mu, hmc_gaugefield* gaugefield){
  set_zero_spinor(b);

  int color = spinor_color(i);
  int spin = spinor_spin(i,color);

  b[spinor_field_element(spin,color,spacepos,timepos)].re = sqrt(2.*kappa);

  return HMC_SUCCESS;
}

//!!CP: LZ should update this...
hmc_error create_point_source_eoprec(hmc_eoprec_spinor_field* be,hmc_eoprec_spinor_field* bo,int i,int spacepos,int timepos,hmc_float kappa, hmc_float mu, hmc_float theta,hmc_float chem_pot_re, hmc_float chem_pot_im, hmc_gaugefield* gaugefield){
  
  hmc_spinor_field* source = new hmc_spinor_field[SPINORFIELDSIZE];

  set_zero_spinor(source);
  int color = spinor_color(i);
  int spin = spinor_spin(i,color);

  source[spinor_field_element(spin,color,spacepos,timepos)].re = hmc_one_f;

  hmc_eoprec_spinor_field evensource[EOPREC_SPINORFIELDSIZE];

  convert_to_eoprec(evensource,bo,source);
  convert_to_kappa_format_eoprec(evensource,kappa);
  convert_to_kappa_format_eoprec(bo,kappa);

  hmc_eoprec_spinor_field spintmp[EOPREC_SPINORFIELDSIZE];
  M_inverse_sitediagonal(spintmp, bo, kappa,mu);
  dslash_eoprec(be,spintmp,gaugefield,kappa,theta, chem_pot_re, chem_pot_im, EVEN);

  for(int n=0;n<EOPREC_SPINORFIELDSIZE;n++) {
    be[n].re = evensource[n].re - be[n].re;
    be[n].im = evensource[n].im - be[n].im;
  }

  delete [] source;
  return HMC_SUCCESS;
}

void spinor_apply_bc(hmc_spinor * in, hmc_float theta){
  for(int n = 0; n<SPINORSIZE; n++){
    hmc_float tmp1 = in[n].re;
    hmc_float tmp2 = in[n].im;
    in[n].re = cos(theta)*tmp1 - sin(theta)*tmp2;
    in[n].im = sin(theta)*tmp1 + cos(theta)*tmp2;
  }
  return; 
}

