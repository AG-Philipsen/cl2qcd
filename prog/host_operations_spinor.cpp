#include "host_operations_spinor.h"

hmc_error set_local_zero_spinor(hmc_spinor* inout){
  for(int j=0; j<SPINORSIZE; j++) {
    inout[j].re = 0;
    inout[j].im = 0;
  }
  return HMC_SUCCESS;
}

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
  for(int j=0; j<SPINORSIZE; j++) {
		inout[j].re *=factor;
		inout[j].im *=factor;
	}
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

hmc_error multiply_spinor_i_factor_gamma5(hmc_spinor* in, hmc_spinor* out, hmc_float factor){
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

//this is defined with an "i" inside!!
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
      out[spinor_element(alpha,c)].re = vec_out[c].re;
      out[spinor_element(alpha,c)].im = vec_out[c].im;
    }
  }
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

hmc_error spinors_accumulate(hmc_spinor* inout, hmc_spinor* incr){
  for(int j=0; j<SPINORSIZE; j++) {
    inout[j].re += incr[j].re;
    inout[j].im += incr[j].im;
  }
  return HMC_SUCCESS;
}

//spinout =  (1 + i*gamma_5*mubar)spin_in
void M_diag_local(hmc_spinor* spininout, hmc_float mubar){
	hmc_spinor spintmp[SPINORSIZE];
	multiply_spinor_i_factor_gamma5(spininout,spintmp,mubar);
	spinors_accumulate(spininout,spintmp);
	return;
}


//spinout = U_0*(r-gamma_0)*spinnext + U^dagger_0(x-hat0) * (r+gamma_0)*spinprev
void dslash_0(hmc_spinor* spinnext, hmc_spinor* spinprev, hmc_spinor* spinout, hmc_su3matrix* u, hmc_su3matrix* udagger){
	
// 	hmc_spinor tmp[SPINORSIZE];
// 	multiply_spinor_gamma0(spinnext,tmp);
// 	real_multiply_spinor(tmp,-hmc_one_f);
// 	spinors_accumulate(spinnext,tmp);
// 	su3matrix_times_spinor(u,spinnext,tmp);
// 	spinors_accumulate(spinout,tmp);
// 
// 	multiply_spinor_gamma0(spinprev,tmp);
// 	spinors_accumulate(spinprev,tmp);
// 	su3matrix_times_spinor(udagger,spinprev,tmp);
// 	spinors_accumulate(spinout,tmp);
	
	spinprojectproduct_gamma0(u,spinnext,-hmc_one_f);
	spinors_accumulate(spinout,spinnext);

	spinprojectproduct_gamma0(udagger,spinprev,hmc_one_f);
	spinors_accumulate(spinout,spinprev);
	
	return;
}

// spinout += U_1*(r-gamma_1)*spinnext + U^dagger_1(x-hat1) * (r+gamma_1)*spinprev
void dslash_1(hmc_spinor* spinnext, hmc_spinor* spinprev, hmc_spinor* spinout, hmc_su3matrix* u, hmc_su3matrix* udagger){
	
// 	hmc_spinor tmp[SPINORSIZE];
// 	multiply_spinor_gamma1(spinnext,tmp);
// 	real_multiply_spinor(tmp,-hmc_one_f);
// 	spinors_accumulate(spinnext,tmp);
// 	su3matrix_times_spinor(u,spinnext,tmp);
// 	spinors_accumulate(spinout,tmp);
// 
// 	multiply_spinor_gamma1(spinprev,tmp);
// 	spinors_accumulate(spinprev,tmp);
// 	su3matrix_times_spinor(udagger,spinprev,tmp);
// 	spinors_accumulate(spinout,tmp);
	
		spinprojectproduct_gamma1(u,spinnext,-hmc_one_f);
	spinors_accumulate(spinout,spinnext);

	spinprojectproduct_gamma1(udagger,spinprev,hmc_one_f);
	spinors_accumulate(spinout,spinprev);
	
	return;
}

// spinout += U_2*(r-gamma_2)*spinnext + U^dagger_2(x-hat2) * (r+gamma_2)*spinprev
void dslash_2(hmc_spinor* spinnext, hmc_spinor* spinprev, hmc_spinor* spinout, hmc_su3matrix* u, hmc_su3matrix* udagger){
	
// 	hmc_spinor tmp[SPINORSIZE];
// 	multiply_spinor_gamma2(spinnext,tmp);
// 	real_multiply_spinor(tmp,-hmc_one_f);
// 	spinors_accumulate(spinnext,tmp);
// 	su3matrix_times_spinor(u,spinnext,tmp);
// 	spinors_accumulate(spinout,tmp);
// 
// 	multiply_spinor_gamma2(spinprev,tmp);
// 	spinors_accumulate(spinprev,tmp);
// 	su3matrix_times_spinor(udagger,spinprev,tmp);
// 	spinors_accumulate(spinout,tmp);
	
		spinprojectproduct_gamma2(u,spinnext,-hmc_one_f);
	spinors_accumulate(spinout,spinnext);

	spinprojectproduct_gamma2(udagger,spinprev,hmc_one_f);
	spinors_accumulate(spinout,spinprev);
	
	return;
}

// spinout += U_3*(r-gamma_3)*spinnext + U^dagger_3(x-hat3) * (r+gamma_3)*spinprev
void dslash_3(hmc_spinor* spinnext, hmc_spinor* spinprev, hmc_spinor* spinout, hmc_su3matrix* u, hmc_su3matrix* udagger){
	
// 	hmc_spinor tmp[SPINORSIZE];
// 	multiply_spinor_gamma3(spinnext,tmp);
// 	real_multiply_spinor(tmp,-hmc_one_f);
// 	spinors_accumulate(spinnext,tmp);
// 	su3matrix_times_spinor(u,spinnext,tmp);
// 	spinors_accumulate(spinout,tmp);
// 
// 	multiply_spinor_gamma3(spinprev,tmp);
// 	spinors_accumulate(spinprev,tmp);
// 	su3matrix_times_spinor(udagger,spinprev,tmp);
// 	spinors_accumulate(spinout,tmp);
	
	spinprojectproduct_gamma3(u,spinnext,-hmc_one_f);
	spinors_accumulate(spinout,spinnext);

	spinprojectproduct_gamma3(udagger,spinprev,hmc_one_f);
	spinors_accumulate(spinout,spinprev);
	
	return;
}

void gamma_5_spinor(hmc_full_spinor inout){
	for(int i = 0; i<NC; i++){
		for(int j = 2; j<NDIM; j++){
			inout[spinor_element(j, i)].re *= -1.;
			inout[spinor_element(j, i)].im *= -1.;
		}}
	
	/** @todo CP: this was originally intended to be used, but can only be applied in this way if the spinorstructure is changed!!*/
// 	su3_vector_times_minusone(&(inout[2*NC]));
// 	su3_vector_times_minusone(&(inout[3*NC]));
}

void su3_vector_times_minusone(hmc_su3vector inout){
	(inout[0]).re *= -1.;
	(inout[0]).im *= -1.;
	(inout[1]).re *= -1.;
	(inout[1]).im *= -1.;
	(inout[2]).re *= -1.;
	(inout[2]).im *= -1.;
}

void su3_vector_acc(hmc_su3vector in, hmc_su3vector out){
	(out[0]).re += (in[0]).re;
	(out[0]).im += (in[0]).im;
	(out[1]).re += (in[1]).re;
	(out[1]).im += (in[1]).im;
	(out[2]).re += (in[2]).re;
	(out[2]).im += (in[2]).im;
}

void su3_vector_multiple_add(hmc_su3vector in1, hmc_su3vector in2, hmc_su3vector out){
	(out[0]).re = (in1[0]).re + (in2[0]).re;
	(out[1]).re = (in1[1]).re + (in2[1]).re;
	(out[2]).re = (in1[2]).re + (in2[2]).re;
	(out[0]).im = (in1[0]).im + (in2[0]).im;
	(out[1]).im = (in1[1]).im + (in2[1]).im;
	(out[2]).im = (in1[2]).im + (in2[2]).im;
}

void spinproj_gamma1_a(hmc_full_spinor spin, hmc_su3vector out, hmc_float sign){
	for(int c = 0; c<NC; c++){
		out[c].re = spin[spinor_element(0,c)].re - sign*spin[spinor_element(3,c)].im;
  	out[c].im = spin[spinor_element(0,c)].im + sign*spin[spinor_element(3,c)].re;
  }
}

void spinproj_gamma1_b(hmc_full_spinor spin, hmc_su3vector out, hmc_float sign){
	for(int c = 0; c<NC; c++){
		out[c].re = spin[spinor_element(1,c)].re - sign*spin[spinor_element(2,c)].im;
  	out[c].im = spin[spinor_element(1,c)].im + sign*spin[spinor_element(2,c)].re;
  }
}

void spinproj_gamma2_a(hmc_full_spinor spin, hmc_su3vector out, hmc_float sign){
	for(int c = 0; c<NC; c++){
		out[c].re = spin[spinor_element(0,c)].re - sign*spin[spinor_element(3,c)].re;
  	out[c].im = spin[spinor_element(0,c)].im - sign*spin[spinor_element(3,c)].im;
  }
}

void spinproj_gamma2_b(hmc_full_spinor spin, hmc_su3vector out, hmc_float sign){
	for(int c = 0; c<NC; c++){
		out[c].re = spin[spinor_element(1,c)].re + sign*spin[spinor_element(2,c)].re;
  	out[c].im = spin[spinor_element(1,c)].im + sign*spin[spinor_element(2,c)].im;
  }
}

void spinproj_gamma3_a(hmc_full_spinor spin, hmc_su3vector out, hmc_float sign){
	for(int c = 0; c<NC; c++){
		out[c].re = spin[spinor_element(0,c)].re - sign*spin[spinor_element(2,c)].im;
  	out[c].im = spin[spinor_element(0,c)].im + sign*spin[spinor_element(2,c)].re;
  }
}

void spinproj_gamma3_b(hmc_full_spinor spin, hmc_su3vector out, hmc_float sign){
	for(int c = 0; c<NC; c++){
		out[c].re = spin[spinor_element(1,c)].re + sign*spin[spinor_element(3,c)].im;
  	out[c].im = spin[spinor_element(1,c)].im - sign*spin[spinor_element(3,c)].re;
  }
}


void spinproj_gamma0_a(hmc_full_spinor spin, hmc_su3vector out, hmc_float sign){
	for(int c = 0; c<NC; c++){
		out[c].re = spin[spinor_element(0,c)].re + sign*spin[spinor_element(2,c)].re;
  	out[c].im = spin[spinor_element(0,c)].im + sign*spin[spinor_element(2,c)].im;
  }
}

void spinproj_gamma0_b(hmc_full_spinor spin, hmc_su3vector out, hmc_float sign){
	for(int c = 0; c<NC; c++){
		out[c].re = spin[spinor_element(1,c)].re + sign*spin[spinor_element(3,c)].re;
  	out[c].im = spin[spinor_element(1,c)].im + sign*spin[spinor_element(3,c)].im;
  }
}


////////////////////////////////////////////////////////////////////
//our original gamma_matrices:

/*

void spinproj_gamma0_a(hmc_full_spinor spin, hmc_su3vector out, hmc_float sign){
	for(int c = 0; c<NC; c++){
		out[c].re = spin[spinor_element(0,c)].re - sign*spin[spinor_element(3,c)].im;
  	out[c].im = spin[spinor_element(0,c)].im + sign*spin[spinor_element(3,c)].re;
  }
}

void spinproj_gamma0_b(hmc_full_spinor spin, hmc_su3vector out, hmc_float sign){
	for(int c = 0; c<NC; c++){
		out[c].re = spin[spinor_element(1,c)].re - sign*spin[spinor_element(2,c)].im;
  	out[c].im = spin[spinor_element(1,c)].im + sign*spin[spinor_element(2,c)].re;
  }
}

void spinproj_gamma1_a(hmc_full_spinor spin, hmc_su3vector out, hmc_float sign){
	for(int c = 0; c<NC; c++){
		out[c].re = spin[spinor_element(0,c)].re - sign*spin[spinor_element(3,c)].re;
  	out[c].im = spin[spinor_element(0,c)].im - sign*spin[spinor_element(3,c)].im;
  }
}

void spinproj_gamma1_b(hmc_full_spinor spin, hmc_su3vector out, hmc_float sign){
	for(int c = 0; c<NC; c++){
		out[c].re = spin[spinor_element(1,c)].re + sign*spin[spinor_element(2,c)].re;
  	out[c].im = spin[spinor_element(1,c)].im + sign*spin[spinor_element(2,c)].im;
  }
}

void spinproj_gamma2_a(hmc_full_spinor spin, hmc_su3vector out, hmc_float sign){
	for(int c = 0; c<NC; c++){
		out[c].re = spin[spinor_element(0,c)].re - sign*spin[spinor_element(2,c)].im;
  	out[c].im = spin[spinor_element(0,c)].im + sign*spin[spinor_element(2,c)].re;
  }
}

void spinproj_gamma2_b(hmc_full_spinor spin, hmc_su3vector out, hmc_float sign){
	for(int c = 0; c<NC; c++){
		out[c].re = spin[spinor_element(1,c)].re + sign*spin[spinor_element(3,c)].im;
  	out[c].im = spin[spinor_element(1,c)].im - sign*spin[spinor_element(3,c)].re;
  }
}


void spinproj_gamma3_a(hmc_full_spinor spin, hmc_su3vector out, hmc_float sign){
	for(int c = 0; c<NC; c++){
		out[c].re = spin[spinor_element(0,c)].re + sign*spin[spinor_element(2,c)].re;
  	out[c].im = spin[spinor_element(0,c)].im + sign*spin[spinor_element(2,c)].im;
  }
}

void spinproj_gamma3_b(hmc_full_spinor spin, hmc_su3vector out, hmc_float sign){
	for(int c = 0; c<NC; c++){
		out[c].re = spin[spinor_element(1,c)].re + sign*spin[spinor_element(3,c)].re;
  	out[c].im = spin[spinor_element(1,c)].im + sign*spin[spinor_element(3,c)].im;
  }
}

*/

//////////////////////////////////////////////////////////////////////


//calculates the trace of i times generator times 3x3-matrix and stores this in a su3-algebraelement
//now using structs
void tr_lambda_u(hmc_3x3matrix in, hmc_algebraelement2 out){
				(out).e0 = ( -(in)[1][0].im - (in)[0][1].im);
				(out).e1 = (+(in)[1][0].re-(in)[0][1].re);
				(out).e2 = (-(in)[0][0].im+(in)[1][1].im);
				(out).e3 = (-(in)[2][0].im-(in)[0][2].im);
				(out).e4 = (+(in)[2][0].re-(in)[0][2].re);
				(out).e5 = (-(in)[2][1].im-(in)[1][2].im);
				(out).e6 = (+(in)[2][1].re-(in)[1][2].re);
				(out).e7 = (-(in)[0][0].im-(in)[1][1].im + 2.0*(in)[2][2].im)*0.577350269189625;
	
}

//deprecated version without struct
///////////////////////////////////////////////////////////////
/*
void tr_lambda_u(hmc_3x3matrix in, hmc_algebraelement out){
				(out)[0] = ( -(in)[1][0].im - (in)[0][1].im);
				(out)[1] = (+(in)[1][0].re-(in)[0][1].re);
				(out)[2] = (-(in)[0][0].im+(in)[1][1].im);
				(out)[3] = (-(in)[2][0].im-(in)[0][2].im);
				(out)[4] = (+(in)[2][0].re-(in)[0][2].re);
				(out)[5] = (-(in)[2][1].im-(in)[1][2].im);
				(out)[6] = (+(in)[2][1].re-(in)[1][2].re);
				(out)[7] = (-(in)[0][0].im-(in)[1][1].im + 2.0*(in)[2][2].im)*0.577350269189625;
}
*/


//calculates the Dirac-Trace of the matrix resulting from multiplying v*u^dagger + w*x^dagger, where u, v, w, x are SU(3)-vectors
//	the result is a 3x3-matrix
void tr_v_times_u_dagger(hmc_su3vector v, hmc_su3vector u, hmc_su3vector w, hmc_su3vector x, hmc_3x3matrix out){
	for(int a = 0; a<3; a++){
		for(int b = 0; b<3; b++){
			(out[a][b]).re = (v[a]).re*(u[b]).re + (v[a]).im*(u[b]).im + (w[a]).re*(x[b]).re + (w[a]).im*(x[b]).im;
			(out[a][b]).im = (v[a]).im*(u[b]).re - (v[a]).re*(u[b]).im + (w[a]).im*(x[b]).re - (w[a]).re*(x[b]).im;;
		}}
	/*
	 #define _vector_tensor_vector_add(t, u, v, w, z) \
   (t).c00.re=(u).c0.re*(v).c0.re+(u).c0.im*(v).c0.im + (w).c0.re*(z).c0.re+(w).c0.im*(z).c0.im; \
   (t).c00.im=(u).c0.im*(v).c0.re-(u).c0.re*(v).c0.im + (w).c0.im*(z).c0.re-(w).c0.re*(z).c0.im; \
   (t).c01.re=(u).c0.re*(v).c1.re+(u).c0.im*(v).c1.im + (w).c0.re*(z).c1.re+(w).c0.im*(z).c1.im; \
   (t).c01.im=(u).c0.im*(v).c1.re-(u).c0.re*(v).c1.im + (w).c0.im*(z).c1.re-(w).c0.re*(z).c1.im; \
   (t).c02.re=(u).c0.re*(v).c2.re+(u).c0.im*(v).c2.im + (w).c0.re*(z).c2.re+(w).c0.im*(z).c2.im; \
   (t).c02.im=(u).c0.im*(v).c2.re-(u).c0.re*(v).c2.im + (w).c0.im*(z).c2.re-(w).c0.re*(z).c2.im; \
   (t).c10.re=(u).c1.re*(v).c0.re+(u).c1.im*(v).c0.im + (w).c1.re*(z).c0.re+(w).c1.im*(z).c0.im; \
   (t).c10.im=(u).c1.im*(v).c0.re-(u).c1.re*(v).c0.im + (w).c1.im*(z).c0.re-(w).c1.re*(z).c0.im; \
   (t).c11.re=(u).c1.re*(v).c1.re+(u).c1.im*(v).c1.im + (w).c1.re*(z).c1.re+(w).c1.im*(z).c1.im; \
   (t).c11.im=(u).c1.im*(v).c1.re-(u).c1.re*(v).c1.im + (w).c1.im*(z).c1.re-(w).c1.re*(z).c1.im; \
   (t).c12.re=(u).c1.re*(v).c2.re+(u).c1.im*(v).c2.im + (w).c1.re*(z).c2.re+(w).c1.im*(z).c2.im; \
   (t).c12.im=(u).c1.im*(v).c2.re-(u).c1.re*(v).c2.im + (w).c1.im*(z).c2.re-(w).c1.re*(z).c2.im; \
   (t).c20.re=(u).c2.re*(v).c0.re+(u).c2.im*(v).c0.im + (w).c2.re*(z).c0.re+(w).c2.im*(z).c0.im; \
   (t).c20.im=(u).c2.im*(v).c0.re-(u).c2.re*(v).c0.im + (w).c2.im*(z).c0.re-(w).c2.re*(z).c0.im; \
   (t).c21.re=(u).c2.re*(v).c1.re+(u).c2.im*(v).c1.im + (w).c2.re*(z).c1.re+(w).c2.im*(z).c1.im; \
   (t).c21.im=(u).c2.im*(v).c1.re-(u).c2.re*(v).c1.im + (w).c2.im*(z).c1.re-(w).c2.re*(z).c1.im; \
   (t).c22.re=(u).c2.re*(v).c2.re+(u).c2.im*(v).c2.im + (w).c2.re*(z).c2.re+(w).c2.im*(z).c2.im; \
   (t).c22.im=(u).c2.im*(v).c2.re-(u).c2.re*(v).c2.im + (w).c2.im*(z).c2.re-(w).c2.re*(z).c2.im; 
	 */
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////
// deprecated functions
///////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
//spinout =  (3 + 4*kappa^2*mu^2)spin_in
void MdaggerM_diag_local(hmc_spinor* spininout, hmc_float kappa, hmc_float mu){
	hmc_float twistfactor = 1. + 8.*kappa*kappa + 4.*kappa*kappa*mu*mu;
	real_multiply_spinor(spininout,twistfactor);
	return;
}


//spinout = U_0*(r+gamma_0)*spinnext + U^dagger_0(x-hat0) * (r-gamma_0)*spinprev
void ddaggerslash_0(hmc_spinor* spinnext, hmc_spinor* spinprev, hmc_spinor* spinout, hmc_su3matrix* u, hmc_su3matrix* udagger){
	
// 	hmc_spinor tmp[SPINORSIZE];
// 	multiply_spinor_gamma0(spinnext,tmp);
// 	spinors_accumulate(spinnext,tmp);
// 	su3matrix_times_spinor(u,spinnext,tmp);
// 	spinors_accumulate(spinout,tmp);
// 
// 	multiply_spinor_gamma0(spinprev,tmp);
// 	real_multiply_spinor(tmp,-hmc_one_f);
// 	spinors_accumulate(spinprev,tmp);
// 	su3matrix_times_spinor(udagger,spinprev,tmp);
// 	spinors_accumulate(spinout,tmp);
	
	spinprojectproduct_gamma0(u,spinnext,hmc_one_f);
	spinors_accumulate(spinout,spinnext);

	spinprojectproduct_gamma0(udagger,spinprev,-hmc_one_f);
	spinors_accumulate(spinout,spinprev);
	
	return;
}

// spinout += U_1*(r+gamma_1)*spinnext + U^dagger_1(x-hat1) * (r-gamma_1)*spinprev
void ddaggerslash_1(hmc_spinor* spinnext, hmc_spinor* spinprev, hmc_spinor* spinout, hmc_su3matrix* u, hmc_su3matrix* udagger){
	
// 	hmc_spinor tmp[SPINORSIZE];
// 	multiply_spinor_gamma1(spinnext,tmp);
// 	spinors_accumulate(spinnext,tmp);
// 	su3matrix_times_spinor(u,spinnext,tmp);
// 	spinors_accumulate(spinout,tmp);
// 
// 	multiply_spinor_gamma1(spinprev,tmp);
// 	real_multiply_spinor(tmp,-hmc_one_f);
// 	spinors_accumulate(spinprev,tmp);
// 	su3matrix_times_spinor(udagger,spinprev,tmp);
// 	spinors_accumulate(spinout,tmp);
	
	spinprojectproduct_gamma1(u,spinnext,hmc_one_f);
	spinors_accumulate(spinout,spinnext);

	spinprojectproduct_gamma1(udagger,spinprev,-hmc_one_f);
	spinors_accumulate(spinout,spinprev);
	
	return;
}

// spinout += U_2*(r+gamma_2)*spinnext + U^dagger_2(x-hat2) * (r-gamma_2)*spinprev
void ddaggerslash_2(hmc_spinor* spinnext, hmc_spinor* spinprev, hmc_spinor* spinout, hmc_su3matrix* u, hmc_su3matrix* udagger){
	
// 	hmc_spinor tmp[SPINORSIZE];
// 	multiply_spinor_gamma2(spinnext,tmp);
// 	spinors_accumulate(spinnext,tmp);
// 	su3matrix_times_spinor(u,spinnext,tmp);
// 	spinors_accumulate(spinout,tmp);
// 
// 	multiply_spinor_gamma2(spinprev,tmp);
// 	real_multiply_spinor(tmp,-hmc_one_f);
// 	spinors_accumulate(spinprev,tmp);
// 	su3matrix_times_spinor(udagger,spinprev,tmp);
// 	spinors_accumulate(spinout,tmp);

	spinprojectproduct_gamma2(u,spinnext,hmc_one_f);
	spinors_accumulate(spinout,spinnext);

	spinprojectproduct_gamma2(udagger,spinprev,-hmc_one_f);
	spinors_accumulate(spinout,spinprev);
	return;
}

// spinout += U_3*(r+gamma_3)*spinnext + U^dagger_3(x-hat3) * (r-gamma_3)*spinprev
void ddaggerslash_3(hmc_spinor* spinnext, hmc_spinor* spinprev, hmc_spinor* spinout, hmc_su3matrix* u, hmc_su3matrix* udagger){
	
// 	hmc_spinor tmp[SPINORSIZE];
// 	multiply_spinor_gamma3(spinnext,tmp);
// 	spinors_accumulate(spinnext,tmp);
// 	su3matrix_times_spinor(u,spinnext,tmp);
// 	spinors_accumulate(spinout,tmp);
// 
// 	multiply_spinor_gamma3(spinprev,tmp);
// 	real_multiply_spinor(tmp,-hmc_one_f);
// 	spinors_accumulate(spinprev,tmp);
// 	su3matrix_times_spinor(udagger,spinprev,tmp);
// 	spinors_accumulate(spinout,tmp);

	spinprojectproduct_gamma3(u,spinnext,hmc_one_f);
	spinors_accumulate(spinout,spinnext);

	spinprojectproduct_gamma3(udagger,spinprev,-hmc_one_f);
	spinors_accumulate(spinout,spinprev);
	
	return;
}

*/
