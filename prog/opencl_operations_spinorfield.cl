//opencl_operations_spinorfield

//eoprec operations
void convert_to_eoprec(hmc_eoprec_spinor_field* even, hmc_eoprec_spinor_field* odd, hmc_spinor_field* in){
  int spacepos, timepos;
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
  return;
}

void convert_from_eoprec(hmc_eoprec_spinor_field* even, hmc_eoprec_spinor_field* odd, hmc_spinor_field* out){
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
  return;
}

void convert_to_kappa_format(__global hmc_spinor_field* inout,hmc_float kappa){
  for(int n=0; n<SPINORFIELDSIZE; n++) {
    inout[n].re *= sqrt(2.*kappa);
    inout[n].im *= sqrt(2.*kappa);
  }
  return;
}

void convert_to_kappa_format_eoprec(__global hmc_eoprec_spinor_field* inout,hmc_float kappa){
  for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
    inout[n].re *= sqrt(2.*kappa);
    inout[n].im *= sqrt(2.*kappa);
  }
  return;
}

void convert_from_kappa_format(__global hmc_spinor_field* in, __global hmc_spinor_field * out,hmc_float kappa){
  for(int n=0; n<SPINORFIELDSIZE; n++) {
    out[n].re = (in[n].re)/sqrt(2.*kappa);
    out[n].im = (in[n].im)/sqrt(2.*kappa);
  }
  return;
}

void convert_from_kappa_format_eoprec(__global hmc_eoprec_spinor_field* in, __global hmc_eoprec_spinor_field * out, hmc_float kappa){
  for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
    out[n].re = (in[n].re)/sqrt(2.*kappa);
    out[n].im = (in[n].im)/sqrt(2.*kappa);
  }
  return;
}


void get_spinor_from_eoprec_field(__global hmc_eoprec_spinor_field* in, hmc_spinor* out, int n_eoprec){
  for(int alpha=0; alpha<NSPIN; alpha++) {
    for(int color=0; color<NC; color++) {
      out[spinor_element(alpha,color)] = in[eoprec_spinor_field_element(alpha,color,n_eoprec)];
    }
  }
  return;
}

void put_spinor_to_eoprec_field(hmc_spinor* in,__global hmc_eoprec_spinor_field* out, int n_eoprec){
  for(int alpha=0; alpha<NSPIN; alpha++) {
    for(int color=0; color<NC; color++) {
      out[eoprec_spinor_field_element(alpha,color,n_eoprec)]=in[spinor_element(alpha,color)];
    }
  }
  return;
}

void get_spinor_from_field(__global hmc_spinor_field* in, hmc_spinor* out, int n, int t){
  for(int alpha=0; alpha<NSPIN; alpha++) {
    for(int color=0; color<NC; color++) {
      out[spinor_element(alpha,color)] = in[spinor_field_element(alpha,color,n,t)];
    }
  }
  return;
}

void put_spinor_to_field(hmc_spinor* in, __global hmc_spinor_field* out, int n, int t){
  for(int alpha=0; alpha<NSPIN; alpha++) {
    for(int color=0; color<NC; color++) {
      out[spinor_field_element(alpha,color,n,t)]=in[spinor_element(alpha,color)];
    }
  }
  return;
}


void create_point_source(hmc_spinor_field* b, int i, int spacepos, int timepos, hmc_float kappa, hmc_float mu, __global hmc_ocl_gaugefield * gaugefield){
  set_zero_spinor(b);

  int color = spinor_color(i);
  int spin = spinor_spin(i,color);

  b[spinor_field_element(spin,color,spacepos,timepos)].re = sqrt(2.*kappa);

  return;
}

//!!CP: LZ should update this...
void create_point_source_eoprec(hmc_eoprec_spinor_field* be,hmc_eoprec_spinor_field* bo,int i,int spacepos,int timepos,hmc_float kappa, hmc_float mu, hmc_float theta,hmc_float chem_pot_re, hmc_float chem_pot_im, __global hmc_ocl_gaugefield* gaugefield){
  
  
  return;
}
