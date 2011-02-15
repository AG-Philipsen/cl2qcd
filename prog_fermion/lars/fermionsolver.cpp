#include "fermionsolver.h"

hmc_error solver(hmc_spinor_field* in, hmc_spinor_field* out, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, int cgmax){

  hmc_eoprec_spinor_field even[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field odd[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field newodd[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field neweven[EOPREC_SPINORFIELDSIZE];

  hmc_eoprec_spinor_field be[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field bo[EOPREC_SPINORFIELDSIZE];

  convert_to_eoprec(even,odd,in);

  convert_to_kappa_format(even,kappa);
  convert_to_kappa_format(odd,kappa);

  create_point_source(be,bo,1,0,0,kappa,mu,gaugefield);

 solve_cg(neweven, even, be, gaugefield, kappa, mu, cgmax);
 //  solve_bicgstab(neweven, even, be, gaugefield, kappa, mu, cgmax);

  hmc_eoprec_spinor_field spintmp[EOPREC_SPINORFIELDSIZE];
  dslash(newodd,neweven,gaugefield,kappa,ODD); //use newood as tmp variable
  inverse_sitediagonal(spintmp,newodd,kappa,mu);
  inverse_sitediagonal(newodd,bo,kappa,mu);
  for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
    newodd[n].re -= spintmp[n].re;
    newodd[n].im -= spintmp[n].im;
  }

  convert_from_kappa_format(neweven,kappa);
  convert_from_kappa_format(newodd,kappa);

  convert_from_eoprec(neweven,newodd,out);

  return HMC_SUCCESS;
}


hmc_error solve_cg(hmc_eoprec_spinor_field* out,hmc_eoprec_spinor_field* in,hmc_eoprec_spinor_field* source,hmc_gaugefield* gaugefield,hmc_float kappa,hmc_float mu,int cgmax){

  //CG according to DeGrand/DeTar, pp.176


  //whats wrong?

  hmc_eoprec_spinor_field rn[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field rnn[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field pn[EOPREC_SPINORFIELDSIZE];

  //init

  for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
    out[n].re = in[n].re;
    out[n].im = in[n].im;
  }

  //iterations:
  for(int iter=0; iter<cgmax; iter++) {

    if(iter%10==0) {
      //fresh start 
      Aee(out,rn,gaugefield,kappa,mu);
      for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
	rn[n].re = source[n].re - rn[n].re;
	rn[n].im = source[n].im - rn[n].im;
      }
      
      printf("true residue squared: %e\n",global_squarenorm_eoprec(rn));
      
      for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
	pn[n].re = rn[n].re;
	pn[n].im = rn[n].im;
      }
    }      

    hmc_eoprec_spinor_field Apn[EOPREC_SPINORFIELDSIZE];
    Aee(pn,Apn,gaugefield,kappa,mu);
    //    hmc_complex num = scalar_product_eoprec(rn,pn);
    //    hmc_complex denom = scalar_product_eoprec(pn,Apn);
    hmc_complex num = scalar_product_eoprec(rn,rn);
    hmc_complex denom = scalar_product_eoprec(pn,Apn);
    hmc_complex alpha = complexdivide(&num,&denom);

    //    printf("alpha: (%f,%f)\n",alpha.re,alpha.im);
  
    for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
      hmc_complex tmp = complexmult(&alpha,&Apn[n]);
      rnn[n].re = rn[n].re - tmp.re;
      rnn[n].im = rn[n].im - tmp.im;
      /*
      if(n==4){
	printf("pn: (%f,%f)\n",pn[n].re,pn[n].im);
	printf("Apn: (%f,%f)\n",Apn[n].re,Apn[n].im);
	printf("alpha*Apn = (%f,%f)\n",tmp.re,tmp.im);
	printf("rn = (%f,%f)\n",rn[n].re,rn[n].im);
	printf("rnn = (%f,%f)\n",rnn[n].re,rnn[n].im);
      }
      */
      tmp = complexmult(&alpha,&pn[n]);
      out[n].re = out[n].re + tmp.re;     
      out[n].im = out[n].im + tmp.im;
    }
    //    printf("new residue: %f\n",global_squarenorm_eoprec(rnn));


    //    hmc_complex tmp1 = scalar_product_eoprec(rnn,rnn);
    //    hmc_complex tmp2 = scalar_product_eoprec(rn,rn);
    //    hmc_complex beta = complexdivide(&tmp1,&tmp2);
    hmc_float beta = global_squarenorm_eoprec(rnn)/global_squarenorm_eoprec(rn);
    //    printf("beta: (%f,%f)\n",beta.re,beta.im);
    for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
      //      hmc_complex tmp = complexmult(&beta,&pn[n]);
      //      pn[n].re = rnn[n].re + tmp.re;
      //      pn[n].im = rnn[n].im + tmp.im;
      pn[n].re = rnn[n].re + beta*pn[n].re;
      pn[n].im = rnn[n].im + beta*pn[n].im;
    }

    hmc_float residue_squared = global_squarenorm_eoprec(rnn);
    hmc_complex check = scalar_product_eoprec(rnn,rnn);
    printf("%d\t%e -- check: %f\n",iter,residue_squared,check.re-residue_squared);

    for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
      rn[n].re = rnn[n].re;
      rn[n].im = rnn[n].im;
    }
        if(iter==120)    return 0;

  
  }

  return HMC_SUCCESS;
}


hmc_error print_spinor(hmc_spinor* in){
  for(int n=0; n<SPINORSIZE; n++) {
    printf("%f ",in[n].re);
  }
  printf("\n");
  for(int n=0; n<SPINORSIZE; n++) {
    printf("%f ",in[n].im);
  }
  printf("\n");
  printf("\n");

  return HMC_SUCCESS;
}

hmc_error solve_bicgstab(hmc_eoprec_spinor_field* out,hmc_eoprec_spinor_field* in,hmc_eoprec_spinor_field* source,hmc_gaugefield* gaugefield,hmc_float kappa,hmc_float mu,int cgmax){

  //BiCGStab according to DeGrand/DeTar, p.180


  // still not correct !!!!!!!!!!!!!!!!!!!!!!


  //init
  for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
    out[n].re = in[n].re;
    out[n].im = in[n].im;
  }


  //iterations:
  for(int iter=0; iter<cgmax; iter++) {
      hmc_eoprec_spinor_field rn[EOPREC_SPINORFIELDSIZE];
      hmc_eoprec_spinor_field rnn[EOPREC_SPINORFIELDSIZE];
      hmc_eoprec_spinor_field pn[EOPREC_SPINORFIELDSIZE];
      hmc_eoprec_spinor_field r0[EOPREC_SPINORFIELDSIZE];

    if(iter%10==0) {
      //fresh start

      //  print_spinor(in);
      Aee(out,r0,gaugefield,kappa,mu);
      //      print_spinor(r0);
      for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
	r0[n].re = source[n].re - r0[n].re;
	r0[n].im = source[n].im - r0[n].im;
      }
   
      for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
	rn[n].re = r0[n].re;
	pn[n].re = r0[n].re;
	rn[n].im = r0[n].im;
	pn[n].im = r0[n].im;
      }

      printf("true residue: %f\n",global_squarenorm_eoprec(r0));

    }

    hmc_eoprec_spinor_field Apn[EOPREC_SPINORFIELDSIZE];
    Aee(pn,Apn,gaugefield,kappa,mu);
    hmc_complex num = scalar_product_eoprec(r0,rn);
    hmc_complex denom = scalar_product_eoprec(r0,Apn);
    hmc_complex alpha = complexdivide(&num,&denom);

    //    printf("alpha: (%f,%f)\n",alpha.re,alpha.im);
  
    hmc_eoprec_spinor_field s[EOPREC_SPINORFIELDSIZE];
    for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
      hmc_complex tmp = complexmult(&alpha,&Apn[n]);
      s[n].re = rn[n].re - tmp.re;
      s[n].im = rn[n].im - tmp.im;
    }

    hmc_eoprec_spinor_field t[EOPREC_SPINORFIELDSIZE];
    Aee(s,t,gaugefield,kappa,mu);
    num = scalar_product_eoprec(t,s);
    denom = scalar_product_eoprec(t,t);
    //   hmc_float rdenom = global_squarenorm_eoprec(t);
    //hmcc_complex omega;
    //    omega.re = num.re/rdenom;
    //    omega.im = num.im/rdenom;
    hmc_complex omega = complexdivide(&num,&denom);
    //    printf("omega: (%f,%f)\n",omega.re,omega.im);
    

    for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
      hmc_complex tmp = complexmult(&omega,&t[n]);
      rnn[n].re = s[n].re - tmp.re;
      rnn[n].im = s[n].im - tmp.im;
      tmp = complexmult(&omega,&s[n]);
      hmc_complex tmp2 = complexmult(&alpha,&pn[n]);
      out[n].re = out[n].re + tmp.re + tmp2.re;
      out[n].im = out[n].im + tmp.im + tmp2.im;
    }
    
    hmc_complex tmp1 = scalar_product_eoprec(r0,rnn);
    hmc_complex tmp2 = scalar_product_eoprec(r0,rn);
    hmc_complex factor1 = complexdivide(&tmp1,&tmp2);
    hmc_complex factor2 = complexdivide(&alpha,&omega);
    hmc_complex beta = complexmult(&factor1,&factor2);
    //   printf("beta: (%f,%f)\n",beta.re,beta.im);

    for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
      hmc_complex omegabeta = complexmult(&omega,&beta);
      hmc_complex tmp1 = complexmult(&omegabeta,&Apn[n]);
      hmc_complex tmp2 = complexmult(&beta,&pn[n]);
      pn[n].re = rnn[n].re + tmp2.re - tmp1.re;
      pn[n].im = rnn[n].im + tmp2.im - tmp1.im;
    }

    hmc_float residue_squared = global_squarenorm_eoprec(rnn);
    printf("%d\t%e\n",iter,residue_squared);

    for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
      rn[n].re = rnn[n].re;
      rn[n].im = rnn[n].im;
    }
  
    if(iter==200) return 0;

  }

  return HMC_SUCCESS;
}


hmc_error create_point_source(hmc_eoprec_spinor_field* be,hmc_eoprec_spinor_field* bo,int i,int spacepos,int timepos,hmc_float kappa, hmc_float mu, hmc_gaugefield* gaugefield){
  hmc_spinor_field source[SPINORFIELDSIZE];
  set_zero_spinor(source);

  int color = spinor_color(i);
  int spin = spinor_spin(i,color);

  source[spinor_field_element(spin,color,spacepos,timepos)].re = hmc_one_f;

  hmc_eoprec_spinor_field evensource[EOPREC_SPINORFIELDSIZE];

  convert_to_eoprec(evensource,bo,source);
  convert_to_kappa_format(evensource,kappa);
  convert_to_kappa_format(bo,kappa);

  hmc_eoprec_spinor_field spintmp[EOPREC_SPINORFIELDSIZE];
  inverse_sitediagonal(spintmp, bo, kappa,mu);
  dslash(be,spintmp,gaugefield,kappa,EVEN);

  for(int n=0;n<EOPREC_SPINORFIELDSIZE;n++) {
    be[n].re = evensource[n].re - be[n].re;
    be[n].im = evensource[n].im - be[n].im;
  }

  return HMC_SUCCESS;
}

/*
hmc_error Aee(hmc_eoprec_spinor_field* in, hmc_eoprec_spinor_field* out, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu){
  for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
    out[n].re = 6.321*in[n].re+in[n].im*0.1;
    out[n].im = -in[n].im;
  }

  out[4].re=out[1].im-.3*out[2].re;
  out[0].im=out[3].re+out[0].re;

  return HMC_SUCCESS;
}
*/


hmc_error Aee(hmc_eoprec_spinor_field* in, hmc_eoprec_spinor_field* out, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu){
  hmc_eoprec_spinor_field spintmp1[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field spintmp2[EOPREC_SPINORFIELDSIZE];

  dslash(spintmp1,in,gaugefield,kappa,ODD); // D_oe
  inverse_sitediagonal(spintmp2,spintmp1,kappa,mu); // R_o^(-1)
  dslash(out,spintmp2,gaugefield,kappa,EVEN); // D_eo
  sitediagonal(spintmp1,in,kappa,mu); //R_e

  for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
    out[n].re = spintmp1[n].re - out[n].re;
    out[n].im = spintmp1[n].im - out[n].im;
  }

  return HMC_SUCCESS;
}



hmc_error dslash(hmc_eoprec_spinor_field* out, hmc_eoprec_spinor_field* in, hmc_gaugefield* gaugefield, hmc_float kappa, int evenodd){
  for(int n=0; n<VOL4D/2; n++) {

    // spinout = U_0*(r-gamma_0)*spinnext + U^dagger_0(x-hat0) * (r+gamma_0)*spinprev
    hmc_spinor spinout[SPINORSIZE];
    set_local_zero_spinor(spinout);    

    int ns = get_nspace_from_eoprecindex(n,evenodd);
    int nt = get_ntime_from_eoprecindex(n,evenodd);
  
    int next = (nt+1)%NTIME;
    int prev = (nt-1+NTIME)%NTIME;
    //implement APBC !!!
    int neo_next = get_n_eoprec(next,ns);
    int neo_prev = get_n_eoprec(prev,ns);

    hmc_spinor spinnext[SPINORSIZE];
    hmc_spinor spinprev[SPINORSIZE];
    get_spinor_from_eoprec_field(in,spinnext,neo_next);
    get_spinor_from_eoprec_field(in,spinprev,neo_prev);

    hmc_su3matrix u;
    get_su3matrix(&u,gaugefield,ns,nt,0);
    hmc_su3matrix udagger;
    get_su3matrix(&udagger,gaugefield,ns,prev,0);
    adjoin_su3matrix(&udagger);

    spinprojectproduct_gamma0(&u,spinnext,-hmc_one_f);
    spinors_accumulate(spinout,spinnext);

    spinprojectproduct_gamma0(&udagger,spinprev,hmc_one_f);
    spinors_accumulate(spinout,spinprev);
  
    //    printf("mu=0 %e\n",spinor_squarenorm(spinout));

    // spinout += U_1*(r-gamma_1)*spinnext + U^dagger_1(x-hat1) * (r+gamma_1)*spinprev
    next = get_neighbor(ns,1);
    prev = get_lower_neighbor(ns,1);
    neo_next = get_n_eoprec(nt,next);
    neo_prev = get_n_eoprec(nt,prev);
    get_spinor_from_eoprec_field(in,spinnext,neo_next);
    get_spinor_from_eoprec_field(in,spinprev,neo_prev);
    
    get_su3matrix(&u,gaugefield,ns,nt,1);
    get_su3matrix(&udagger,gaugefield,prev,nt,1);
    adjoin_su3matrix(&udagger);

    spinprojectproduct_gamma1(&u,spinnext,-hmc_one_f);
    spinors_accumulate(spinout,spinnext);

    spinprojectproduct_gamma1(&udagger,spinprev,hmc_one_f);
    spinors_accumulate(spinout,spinprev);

    //    printf("mu=1 %e\n",spinor_squarenorm(spinout));

    // spinout += U_2*(r-gamma_2)*spinnext + U^dagger_2(x-hat2) * (r+gamma_2)*spinprev
    next = get_neighbor(ns,2);
    prev = get_lower_neighbor(ns,2);
    neo_next = get_n_eoprec(nt,next);
    neo_prev = get_n_eoprec(nt,prev);
    get_spinor_from_eoprec_field(in,spinnext,neo_next);
    get_spinor_from_eoprec_field(in,spinprev,neo_prev);
    
    get_su3matrix(&u,gaugefield,ns,nt,2);
    get_su3matrix(&udagger,gaugefield,prev,nt,2);
    adjoin_su3matrix(&udagger);
      
    spinprojectproduct_gamma2(&u,spinnext,-hmc_one_f);
    spinors_accumulate(spinout,spinnext);

    spinprojectproduct_gamma2(&udagger,spinprev,hmc_one_f);
    spinors_accumulate(spinout,spinprev);
    
    //    printf("mu=2 %e\n",spinor_squarenorm(spinout));
    
    // spinout += U_3*(r-gamma_3)*spinnext + U^dagger_3(x-hat3) * (r+gamma_3)*spinprev
    next = get_neighbor(ns,3);
    prev = get_lower_neighbor(ns,3);
    neo_next = get_n_eoprec(nt,next);
    neo_prev = get_n_eoprec(nt,prev);

    get_spinor_from_eoprec_field(in,spinnext,neo_next);
    get_spinor_from_eoprec_field(in,spinprev,neo_prev);
    
    get_su3matrix(&u,gaugefield,ns,nt,3);
    get_su3matrix(&udagger,gaugefield,prev,nt,3);
    adjoin_su3matrix(&udagger);
      
    spinprojectproduct_gamma3(&u,spinnext,-hmc_one_f);
    spinors_accumulate(spinout,spinnext);

    spinprojectproduct_gamma3(&udagger,spinprev,hmc_one_f);
    spinors_accumulate(spinout,spinprev);

    //    printf("mu=3 %e\n",spinor_squarenorm(spinout));

    real_multiply_spinor(spinout,-kappa);
    put_spinor_to_eoprec_field(spinout,out,n);
  
  }

  return HMC_SUCCESS;
}



hmc_error sitediagonal(hmc_eoprec_spinor_field* out, hmc_eoprec_spinor_field* in, hmc_float kappa, hmc_float mu){
  for(int n=0; n<VOL4D/2; n++) {
    hmc_spinor tmp[SPINORSIZE];
    hmc_spinor spinout[SPINORSIZE];
    get_spinor_from_eoprec_field(in,spinout,n);
    hmc_float twistfactor = 2*kappa*mu;
    multiply_spinor_factor_gamma5(spinout,tmp,twistfactor);
    spinors_accumulate(spinout,tmp);
    put_spinor_to_eoprec_field(spinout,out,n);
  }
  return HMC_SUCCESS;
}

hmc_error inverse_sitediagonal(hmc_eoprec_spinor_field* out, hmc_eoprec_spinor_field* in, hmc_float kappa, hmc_float mu){
  for(int n=0; n<VOL4D/2; n++) {
    hmc_spinor tmp[SPINORSIZE];
    hmc_spinor spinout[SPINORSIZE];
    get_spinor_from_eoprec_field(in,spinout,n);
    hmc_float twistfactor = -2*kappa*mu;
    multiply_spinor_factor_gamma5(spinout,tmp,twistfactor);
    spinors_accumulate(spinout,tmp);
    hmc_float denom = 1 + twistfactor*twistfactor;
    real_multiply_spinor(spinout,1./denom);
    put_spinor_to_eoprec_field(spinout,out,n);
  }
  return HMC_SUCCESS;
}
