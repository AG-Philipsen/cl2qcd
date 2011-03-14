#include "host_operations_spinorfield.h"

hmc_error set_zero_spinorfield(hmc_spinor_field *field) {
  for (int n=0; n<SPINORFIELDSIZE; n++) {
	field[n].re=0;
	field[n].im=0;
  }
  return HMC_SUCCESS;
}

hmc_error set_zero_spinorfield_eoprec(hmc_eoprec_spinor_field *field) {
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
  set_zero_spinorfield(b);

  int color = spinor_color(i);
  int spin = spinor_spin(i,color);

  b[spinor_field_element(spin,color,spacepos,timepos)].re = sqrt(2.*kappa);

  return HMC_SUCCESS;
}

//!!CP: LZ should update this...
hmc_error create_point_source_eoprec(hmc_eoprec_spinor_field* be,hmc_eoprec_spinor_field* bo,int i,int spacepos,int timepos,hmc_float kappa, hmc_float mu, hmc_float theta,hmc_float chem_pot_re, hmc_float chem_pot_im, hmc_gaugefield* gaugefield){
  
  hmc_spinor_field* source = new hmc_spinor_field[SPINORFIELDSIZE];

  set_zero_spinorfield(source);
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
