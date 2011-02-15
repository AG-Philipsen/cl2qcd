#include "simple_fermionsolver.h"

hmc_error simple_solver(hmc_spinor_field* in, hmc_spinor_field* out, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, int cgmax){

  //convert to kappa format
  for(int n=0; n<SPINORFIELDSIZE; n++) {
    in[n].re *= sqrt(2.*kappa);
    in[n].im *= sqrt(2.*kappa);
  }

  hmc_spinor_field b[SPINORFIELDSIZE];
  create_point_source(b,1,0,0,kappa,mu,gaugefield);
  
  simple_solve_cg(in, b, gaugefield, kappa, mu, cgmax);

  //convert from kappa format
  for(int n=0; n<SPINORFIELDSIZE; n++) {
    out[n].re = in[n].re/sqrt(2.*kappa);
    out[n].im = in[n].im/sqrt(2.*kappa);
  }

  return HMC_SUCCESS;
}


hmc_error simple_solve_cg(hmc_spinor_field* inout, hmc_eoprec_spinor_field* source, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, int cgmax){

  //CG according to DeGrand/DeTar, pp.176

  hmc_spinor_field pn[SPINORFIELDSIZE];
  hmc_spinor_field rn[SPINORFIELDSIZE];
  hmc_spinor_field rnn[SPINORFIELDSIZE];

  for(int iter=0; iter<cgmax; iter++) {

    if(iter%20==0) {
      //fresh start
      hmc_spinor_field Mx[SPINORFIELDSIZE];
      M(inout,Mx,gaugefield,kappa,mu);
      for(int n=0; n<SPINORFIELDSIZE; n++) {
	rn[n].re = source[n].re - Mx[n].re;
	rn[n].im = source[n].im - Mx[n].im;
      }

      hmc_complex trueres = scalar_product(rn,rn);
      printf("true residue: (%e,%e)\n",trueres.re,trueres.im);

      for(int n=0; n<SPINORFIELDSIZE; n++) {
	pn[n].re = rn[n].re;
	pn[n].im = rn[n].im;
      }
      
    }

    hmc_spinor_field Mpn[SPINORFIELDSIZE];
    M(pn,Mpn,gaugefield,kappa,mu);
    //    hmc_complex check1 = scalar_product(pn,pn);
    //    hmc_complex check2 = scalar_product(Mpn,Mpn);
    //    printf("pn - Mpn :  (%f,%f) - (%f,%f) \n",check1.re,check1.im,check2.re,check2.im);
    hmc_complex tmp1 = scalar_product(rn,rn);
    hmc_complex tmp2 = scalar_product(pn,Mpn);
    hmc_complex alpha = complexdivide(&tmp1,&tmp2);

    //    printf("alpha: (%e,%e)\n",alpha.re,alpha.im);
    
    for(int n=0; n<SPINORFIELDSIZE; n++) {
      hmc_complex mult = complexmult(&alpha,&Mpn[n]);
      rnn[n].re = rn[n].re - mult.re;
      rnn[n].im = rn[n].im - mult.im;
      mult = complexmult(&alpha,&pn[n]);
      inout[n].re = inout[n].re + mult.re;
      inout[n].im = inout[n].im + mult.im;
    }
    
    hmc_complex tmp3 = scalar_product(rnn,rnn);
    hmc_complex tmp4 = scalar_product(rn,rn);
    hmc_float beta = tmp3.re/tmp4.re;
    
    for(int n=0; n<SPINORFIELDSIZE; n++) {
      pn[n].re = rnn[n].re + beta*pn[n].re;
      pn[n].im = rnn[n].im + beta*pn[n].im;
    }
    
    hmc_complex resid = scalar_product(rnn,rnn);
    printf("%d\t%e\n",iter,resid.re);

    for(int n=0; n<SPINORFIELDSIZE; n++) {
      rn[n].re = rnn[n].re;
      rn[n].im = rnn[n].im;
    }

    if((iter+1)%100==0) return 0;

  }

  return HMC_SUCCESS;
}


hmc_error create_point_source(hmc_spinor_field* b, int i, int spacepos, int timepos, hmc_float kappa, hmc_float mu, hmc_gaugefield* gaugefield){
  set_zero_spinor(b);

  int color = spinor_color(i);
  int spin = spinor_spin(i,color);

  b[spinor_field_element(spin,color,spacepos,timepos)].re = sqrt(kappa);

  return HMC_SUCCESS;
}


hmc_error M(hmc_spinor_field* in, hmc_spinor_field* out, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu){
  
  for(int spacepos=0; spacepos<VOLSPACE; spacepos++) {
    for(int timepos=0; timepos<NTIME; timepos++) {
      hmc_spinor spinout[SPINORSIZE];
      hmc_spinor spintmp[SPINORSIZE];
      get_spinor_from_field(in,spinout,spacepos,timepos);
      hmc_float twistfactor = 2*kappa*mu;
      multiply_spinor_factor_gamma5(spinout,spintmp,twistfactor);
      spinors_accumulate(spinout,spintmp);
      put_spinor_to_field(spinout,out,spacepos,timepos);
    }
  }
  
  hmc_spinor_field tmp[SPINORFIELDSIZE];
  dslash(in,tmp,gaugefield,kappa);

  for(int n=0; n<SPINORFIELDSIZE; n++) {
    out[n].re += tmp[n].re;
    out[n].im += tmp[n].im;
  }

  return HMC_SUCCESS;
}



hmc_error dslash(hmc_spinor_field* in, hmc_spinor_field* out, hmc_gaugefield* gaugefield, hmc_float kappa){

  set_local_zero_spinor(out);    

  for(int spacepos=0; spacepos<VOLSPACE; spacepos++) {
    for(int timepos=0; timepos<NTIME; timepos++) {
    // spinout = U_0*(r-gamma_0)*spinnext + U^dagger_0(x-hat0) * (r+gamma_0)*spinprev
      hmc_spinor spinnext[SPINORSIZE];
      hmc_spinor spinprev[SPINORSIZE];
      int next = (timepos+1)%NTIME;
      int prev = (timepos-1+NTIME)%NTIME;
      //implement APBC!!!
      get_spinor_from_field(in, spinnext, spacepos, next);
      get_spinor_from_field(in, spinprev, spacepos, prev);
      hmc_spinor tmp[SPINORSIZE];

      hmc_su3matrix u;
      get_su3matrix(&u,gaugefield,spacepos,timepos,0);
      hmc_su3matrix udagger;
      get_su3matrix(&udagger,gaugefield,spacepos,prev,0);
      adjoin_su3matrix(&udagger);

      multiply_spinor_gamma0(spinnext,tmp);
      real_multiply_spinor(tmp,-hmc_one_f);
      spinors_accumulate(spinnext,tmp);
      su3matrix_times_spinor(&u,spinnext,tmp);
      spinors_accumulate(out,tmp);

      multiply_spinor_gamma0(spinprev,tmp);
      spinors_accumulate(spinprev,tmp);
      su3matrix_times_spinor(&udagger,spinprev,tmp);
      spinors_accumulate(out,tmp);


    // spinout += U_1*(r-gamma_1)*spinnext + U^dagger_1(x-hat1) * (r+gamma_1)*spinprev


    // spinout += U_2*(r-gamma_2)*spinnext + U^dagger_2(x-hat2) * (r+gamma_2)*spinprev
    
    // spinout += U_3*(r-gamma_3)*spinnext + U^dagger_3(x-hat3) * (r+gamma_3)*spinprev

    }
  }


  return HMC_SUCCESS;
}


