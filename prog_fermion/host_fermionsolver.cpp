#include "host_fermionsolver.h"

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

  solve_bicgstab(neweven, even, be, gaugefield, kappa, mu, cgmax);

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

  //init
  for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
    out[n].re = in[n].re;
    out[n].im = in[n].im;
  }

  //BiCGStab according to hep-lat/9404013
  hmc_eoprec_spinor_field* rn = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field* rhat = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field* v = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field* p = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field* s = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field* t = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_complex rho;
  hmc_complex rho_next;
  hmc_complex alpha;
  hmc_complex omega;
  hmc_complex beta;

  hmc_complex tmp1;
  hmc_complex tmp2;
  hmc_complex tmp3;


  for(int iter=0; iter<cgmax; iter++){

    if(iter%iter_refresh==0) {
      //fresh start
      Aee(out,rn,gaugefield,kappa,mu);
      for (int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
	rn[n].re = source[n].re - rn[n].re;
	rn[n].im = source[n].im - rn[n].im;
	rhat[n].re = rn[n].re;
	rhat[n].im = rn[n].im;
      }
      
      printf("true residue squared: %e\n",global_squarenorm_eoprec(rn));

      rho.re = hmc_one_f;
      alpha.re = hmc_one_f;
      omega.re = hmc_one_f;
      rho.im = 0;
      alpha.im = 0;
      omega.im = 0;

      set_zero_eoprec_spinor(v);
      set_zero_eoprec_spinor(p);
    }

    rho_next = scalar_product_eoprec(rhat,rn);
    tmp1 = complexdivide(&rho_next,&rho);
    rho = rho_next;
    tmp2 = complexdivide(&alpha,&omega);
    beta = complexmult(&tmp1,&tmp2);

    for (int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
      tmp1 = complexmult(&omega,&v[n]);
      tmp2.re = p[n].re - tmp1.re;
      tmp2.im = p[n].im - tmp1.im;
      tmp3 = complexmult(&beta,&tmp2);
      p[n].re = rn[n].re + tmp3.re;
      p[n].im = rn[n].im + tmp3.im;
    }

    Aee(p,v,gaugefield,kappa,mu);

    tmp1 = scalar_product_eoprec(rhat,v);
    alpha = complexdivide(&rho,&tmp1);

    for (int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
      tmp1 = complexmult(&alpha,&v[n]);
      s[n].re = rn[n].re - tmp1.re;
      s[n].im = rn[n].im - tmp1.im;
    }

    Aee(s,t,gaugefield,kappa,mu);

    tmp1 = scalar_product_eoprec(t,s);
    tmp2 = scalar_product_eoprec(t,t);

    omega = complexdivide(&tmp1,&tmp2);

    for (int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
      tmp1 = complexmult(&omega,&s[n]);
      tmp2 = complexmult(&alpha,&p[n]);
      tmp3 = complexmult(&omega,&t[n]);
      out[n].re += tmp1.re + tmp2.re;
      out[n].im += tmp1.im + tmp2.im;
      rn[n].re = s[n].re - tmp3.re;
      rn[n].im = s[n].im - tmp3.im;
    }

    hmc_float resid = global_squarenorm_eoprec(rn);
    printf("%d\t%e\n",iter,resid);

    if(resid<epssquare) {
      hmc_eoprec_spinor_field* aux = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
      Aee(out,aux,gaugefield,kappa,mu);
      for (int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
	aux[n].re = source[n].re - aux[n].re;
	aux[n].im = source[n].im - aux[n].im;
      }
      hmc_float trueresid = global_squarenorm_eoprec(aux);
      printf("true residue squared: %e\n",trueresid);

      if(trueresid<epssquare) return HMC_SUCCESS;

      delete [] aux;
    }



  }

  delete [] rn;
  delete [] rhat;
  delete [] v;
  delete [] p;
  delete [] s;
  delete [] t;
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
    out[n].re = 4000*in[n].re + 1000*in[(n+1)%EOPREC_SPINORFIELDSIZE].im - 4000*in[(n-1)%EOPREC_SPINORFIELDSIZE].re +214*in[(n+1)%EOPREC_SPINORFIELDSIZE].re;
    out[n].im = in[n].im;
  }

  //  out[4].re=out[1].im-.3*out[2].re;
  //  out[0].im=out[3].re+out[0].re;

  return HMC_SUCCESS;
}

*/



hmc_error Aee(hmc_eoprec_spinor_field* in, hmc_eoprec_spinor_field* out, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu){
  hmc_eoprec_spinor_field* spintmp1 = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field* spintmp2 = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];

  dslash(spintmp1,in,gaugefield,kappa,ODD); // D_oe
  inverse_sitediagonal(spintmp2,spintmp1,kappa,mu); // R_o^(-1)
  dslash(out,spintmp2,gaugefield,kappa,EVEN); // D_eo
  sitediagonal(spintmp1,in,kappa,mu); //R_e

  for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
    out[n].re = spintmp1[n].re - out[n].re;
    out[n].im = spintmp1[n].im - out[n].im;
  }

  delete [] spintmp1;
  delete [] spintmp2;

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
