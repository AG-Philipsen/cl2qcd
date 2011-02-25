//!!CP: LZ should update this...

#include "host_solver.h"

hmc_error solver(hmc_spinor_field* in, hmc_spinor_field* out, hmc_spinor_field* b, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im, int cgmax){
  //CP: in the end this should be done with a compiler option
  if(use_eo){
    hmc_eoprec_spinor_field even[EOPREC_SPINORFIELDSIZE];
    hmc_eoprec_spinor_field odd[EOPREC_SPINORFIELDSIZE];
    hmc_eoprec_spinor_field newodd[EOPREC_SPINORFIELDSIZE];
    hmc_eoprec_spinor_field neweven[EOPREC_SPINORFIELDSIZE];

    hmc_eoprec_spinor_field be[EOPREC_SPINORFIELDSIZE];
    hmc_eoprec_spinor_field bo[EOPREC_SPINORFIELDSIZE];
    convert_to_eoprec(even,odd,in);
    convert_to_kappa_format_eoprec(even,kappa);
    convert_to_kappa_format_eoprec(odd,kappa);
 
    create_point_source_eoprec(be,bo,1,0,0,kappa,mu,theta, chem_pot_re, chem_pot_im, gaugefield);
    bicgstab_eoprec(neweven, even, be, gaugefield, kappa, mu, theta, chem_pot_re, chem_pot_im, cgmax);

    hmc_eoprec_spinor_field spintmp[EOPREC_SPINORFIELDSIZE];
    dslash_eoprec(newodd,neweven,gaugefield,kappa, theta, chem_pot_re, chem_pot_im, ODD); //use newood as tmp variable
    M_inverse_sitediagonal(spintmp,newodd,kappa,mu);
    M_inverse_sitediagonal(newodd,bo,kappa,mu);
    for(int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
      newodd[n].re -= spintmp[n].re;
      newodd[n].im -= spintmp[n].im;
    }

    convert_from_kappa_format_eoprec(neweven,neweven, kappa);
    convert_from_kappa_format_eoprec(newodd, newodd, kappa);
    convert_from_eoprec(neweven,newodd,out);
  
    return HMC_SUCCESS;
  }
  else{
    convert_to_kappa_format(in, kappa);
    bicgstab(in, b, gaugefield, kappa, mu, theta, chem_pot_re, chem_pot_im, cgmax);
    convert_from_kappa_format(in, out, kappa);
    
    return HMC_SUCCESS;
  }
  
}

hmc_error bicgstab(hmc_spinor_field* inout, hmc_eoprec_spinor_field* source, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im, int cgmax){

  //BiCGStab according to hep-lat/9404013

  hmc_spinor_field* rn = new hmc_spinor_field[SPINORFIELDSIZE];
  hmc_spinor_field* rhat = new hmc_spinor_field[SPINORFIELDSIZE];
  hmc_spinor_field* v = new hmc_spinor_field[SPINORFIELDSIZE];
  hmc_spinor_field* p = new hmc_spinor_field[SPINORFIELDSIZE];
  hmc_spinor_field* s = new hmc_spinor_field[SPINORFIELDSIZE];
  hmc_spinor_field* t = new hmc_spinor_field[SPINORFIELDSIZE];
  hmc_complex rho;
  hmc_complex rho_next;
  hmc_complex alpha;
  hmc_complex omega;
  hmc_complex beta;

  hmc_complex tmp1;
  hmc_complex tmp2;
  hmc_complex one = hmc_complex_one;
  hmc_complex minusone = hmc_complex_minusone;

  for(int iter=0; iter<cgmax; iter++){

    if(iter%iter_refresh==0) {
      //fresh start
      M(inout,rn,gaugefield,kappa,mu, theta, chem_pot_re, chem_pot_im);
      saxpy(rn, source, &one, rn);
      copy_spinor(rn, rhat);

      alpha = hmc_complex_one;
      omega = hmc_complex_one;
      rho = hmc_complex_one;
      set_zero_spinor(v);
      set_zero_spinor(p);
      
      printf("true residue squared: %e\n",global_squarenorm(rn));
    }

    rho_next = scalar_product(rhat,rn);
    tmp1 = complexdivide(&rho_next,&rho);
    rho = rho_next;
    tmp2 = complexdivide(&alpha,&omega);
    beta = complexmult(&tmp1,&tmp2);

    tmp1 = complexmult(&beta,&omega);
    tmp2 = complexmult(&minusone,&tmp1);
    saxsbypz(p, v, rn, &beta, &tmp2, p);

    M(p,v,gaugefield,kappa,mu, theta, chem_pot_re, chem_pot_im);

    tmp1 = scalar_product(rhat,v);
    alpha = complexdivide(&rho,&tmp1);

    saxpy(v, rn, &alpha, s);

    M(s,t,gaugefield,kappa,mu, theta, chem_pot_re, chem_pot_im);

    tmp1 = scalar_product(t,s);
    tmp2 = scalar_product(t,t);
    omega = complexdivide(&(tmp1),&(tmp2));

    saxpy(t, s, &omega, rn);

    saxsbypz(p, s, inout, &alpha, &omega, inout);

    hmc_float resid = global_squarenorm(rn);

    if(resid<epssquare) {
      hmc_spinor_field* aux = new hmc_spinor_field[SPINORFIELDSIZE];
      M(inout,aux,gaugefield,kappa,mu, theta, chem_pot_re, chem_pot_im);
      saxpy(aux, source, &one, aux);
      hmc_float trueresid = global_squarenorm(aux);
      printf("%d\t%e\t%e\n",iter,resid,trueresid);

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


hmc_error bicgstab_eoprec(hmc_eoprec_spinor_field* out,hmc_eoprec_spinor_field* in,hmc_eoprec_spinor_field* source,hmc_gaugefield* gaugefield,hmc_float kappa,hmc_float mu, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im,int cgmax){

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
  hmc_complex one = hmc_complex_one;
  hmc_complex minusone = hmc_complex_minusone;

  //CP: why is this not needed in bicgstab??
  //init
  copy_spinor_eoprec(in, out);

  for(int iter=0; iter<cgmax; iter++){

    if(iter%iter_refresh==0) {
      //fresh start
      Aee(out,rn,gaugefield,kappa,mu, theta, chem_pot_re, chem_pot_im);
      saxpy_eoprec(rn, source, &one, rn);
      copy_spinor_eoprec(rn, rhat);
      
      alpha = hmc_complex_one;
      omega = hmc_complex_one;
      rho = hmc_complex_one;      
      set_zero_eoprec_spinor(v);
      set_zero_eoprec_spinor(p);
      
      printf("true residue squared: %e\n",global_squarenorm_eoprec(rn));
    }

    rho_next = scalar_product_eoprec(rhat,rn);
    tmp1 = complexdivide(&rho_next,&rho);
    rho = rho_next;
    tmp2 = complexdivide(&alpha,&omega);
    beta = complexmult(&tmp1,&tmp2);

    tmp1 = complexmult(&beta,&omega);
    tmp2 = complexmult(&minusone,&tmp1);
    saxsbypz_eoprec(p, v, rn, &beta, &tmp2, p);
    
    Aee(p,v,gaugefield,kappa,mu, theta, chem_pot_re, chem_pot_im);

    tmp1 = scalar_product_eoprec(rhat,v);
    alpha = complexdivide(&rho,&tmp1);

    saxpy_eoprec(v, rn, &alpha, s);

    Aee(s,t,gaugefield,kappa,mu, theta, chem_pot_re, chem_pot_im);

    tmp1 = scalar_product_eoprec(t,s);
    tmp2 = scalar_product_eoprec(t,t);
    omega = complexdivide(&tmp1,&tmp2);

    saxpy_eoprec(t, s, &omega, rn);

    saxsbypz_eoprec(p, s, out, &alpha, &omega, out);
  
    //CP: is this formula right?? There should be no t i guess...
    //        9. xi = xi−1 + αpi + ωis
//     for (int n=0; n<EOPREC_SPINORFIELDSIZE; n++) {
//       tmp1 = complexmult(&omega,&s[n]);
//       tmp2 = complexmult(&alpha,&p[n]);
//       tmp3 = complexmult(&omega,&t[n]);
//       out[n].re += tmp1.re + tmp2.re;
//       out[n].im += tmp1.im + tmp2.im;
//       rn[n].re = s[n].re - tmp3.re;
//       rn[n].im = s[n].im - tmp3.im;
//     }

    hmc_float resid = global_squarenorm_eoprec(rn);
    printf("%d\t%e\n",iter,resid);

    if(resid<epssquare) {
      hmc_eoprec_spinor_field* aux = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
      Aee(out,aux,gaugefield,kappa,mu, theta, chem_pot_re, chem_pot_im);
      saxpy_eoprec(aux, source, &one, aux);
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


hmc_error M(hmc_spinor_field* in, hmc_spinor_field* out, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im){
  
  M_diag(in, out, kappa, mu);    
  hmc_spinor_field tmp[SPINORFIELDSIZE];
  dslash(in,tmp,gaugefield, theta, chem_pot_re, chem_pot_im);

  hmc_complex kappa_cmplx = {kappa, 0.};
  saxpy(tmp, out, &kappa_cmplx, out);

  return HMC_SUCCESS;
}

hmc_error M_diag(hmc_spinor_field* in, hmc_spinor_field* out, hmc_float kappa, hmc_float mu){
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
  return HMC_SUCCESS; 
}

hmc_error M_sitediagonal(hmc_eoprec_spinor_field* out, hmc_eoprec_spinor_field* in, hmc_float kappa, hmc_float mu){
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

hmc_error M_inverse_sitediagonal(hmc_eoprec_spinor_field* out, hmc_eoprec_spinor_field* in, hmc_float kappa, hmc_float mu){
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

hmc_error dslash(hmc_spinor_field* in, hmc_spinor_field* out, hmc_gaugefield* gaugefield, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im){

  for(int spacepos=0; spacepos<VOLSPACE; spacepos++) {
    for(int timepos=0; timepos<NTIME; timepos++) {

      hmc_spinor spinout[SPINORSIZE];
      hmc_spinor spinnext[SPINORSIZE];
      hmc_spinor spinprev[SPINORSIZE];
      hmc_spinor tmp[SPINORSIZE];
      hmc_su3matrix u;
      hmc_su3matrix udagger;
      int next;
      int prev;

      
      
      
      
      set_local_zero_spinor(spinout);   
      
      
      
      
      
    
      //like in host_geometry
      int coord[NDIM];
      coord[0]=0;
      for(int j=1;j<NDIM;j++) coord[j] = get_spacecoord(spacepos,j);
        
      // spinout = U_0*(r-gamma_0)*spinnext + U^dagger_0(x-hat0) * (r+gamma_0)*spinprev
      next = (timepos+1)%NTIME;
      prev = (timepos-1+NTIME)%NTIME;

      
      
      
      get_spinor_from_field(in, spinnext, spacepos, next);
      get_spinor_from_field(in, spinprev, spacepos, prev);

      if(next == 0) spinor_apply_bc(spinnext, theta);
      else if(prev == NTIME-1) spinor_apply_bc(spinprev, theta);
      
      get_su3matrix(&u,gaugefield,spacepos,timepos,0);
      get_su3matrix(&udagger,gaugefield,spacepos,prev,0);
      adjoin_su3matrix(&udagger);

      //update links with chemical potential, this shall be put into compiler option lateron
      gaugefield_apply_chem_pot(&u, &udagger, chem_pot_re, chem_pot_im);
      
      multiply_spinor_gamma0(spinnext,tmp);
      real_multiply_spinor(tmp,-hmc_one_f);
      spinors_accumulate(spinnext,tmp);
      su3matrix_times_spinor(&u,spinnext,tmp);
      spinors_accumulate(spinout,tmp);

      multiply_spinor_gamma0(spinprev,tmp);
      spinors_accumulate(spinprev,tmp);
      su3matrix_times_spinor(&udagger,spinprev,tmp);
      spinors_accumulate(spinout,tmp);

    
      
      
      
      
      
      
      // spinout += U_1*(r-gamma_1)*spinnext + U^dagger_1(x-hat1) * (r+gamma_1)*spinprev
      next = get_neighbor(spacepos,1);
      prev = get_lower_neighbor(spacepos,1);
      
      
      get_spinor_from_field(in, spinnext, next, timepos);
      get_spinor_from_field(in, spinprev, prev, timepos);
      
      get_su3matrix(&u,gaugefield,spacepos,timepos,1);
      get_su3matrix(&udagger,gaugefield,prev,timepos,1);
      adjoin_su3matrix(&udagger);

      if(coord[1] == NSPACE-1) spinor_apply_bc(spinnext, theta);
      else if(coord[1] == 0) spinor_apply_bc(spinprev, theta);
      
      multiply_spinor_gamma1(spinnext,tmp);
      real_multiply_spinor(tmp,-hmc_one_f);
      spinors_accumulate(spinnext,tmp);
      su3matrix_times_spinor(&u,spinnext,tmp);
      spinors_accumulate(spinout,tmp);

      multiply_spinor_gamma1(spinprev,tmp);
      spinors_accumulate(spinprev,tmp);
      su3matrix_times_spinor(&u,spinprev,tmp);
      spinors_accumulate(spinout,tmp);


      
      
      
      
      
      
      // spinout += U_2*(r-gamma_2)*spinnext + U^dagger_2(x-hat2) * (r+gamma_2)*spinprev
      next = get_neighbor(spacepos,2);
      prev = get_lower_neighbor(spacepos,2);
      
      
      get_spinor_from_field(in, spinnext, next, timepos);
      get_spinor_from_field(in, spinprev, prev, timepos);
      
      get_su3matrix(&u,gaugefield,spacepos,timepos,2);
      get_su3matrix(&udagger,gaugefield,prev,timepos,2);
      adjoin_su3matrix(&udagger);
      
      if(coord[2] == NSPACE-1) spinor_apply_bc(spinnext, theta);
      else if(coord[2] == 0) spinor_apply_bc(spinprev, theta);

      multiply_spinor_gamma2(spinnext,tmp);
      real_multiply_spinor(tmp,-hmc_one_f);
      spinors_accumulate(spinnext,tmp);
      su3matrix_times_spinor(&u,spinnext,tmp);
      spinors_accumulate(spinout,tmp);

      multiply_spinor_gamma2(spinprev,tmp);
      spinors_accumulate(spinprev,tmp);
      su3matrix_times_spinor(&u,spinprev,tmp);
      spinors_accumulate(spinout,tmp);    

      
      
      
      
      
      
      
      
      // spinout += U_3*(r-gamma_3)*spinnext + U^dagger_3(x-hat3) * (r+gamma_3)*spinprev
      next = get_neighbor(spacepos,3);
      prev = get_lower_neighbor(spacepos,3);
      
      
      
      get_spinor_from_field(in, spinnext, next, timepos);
      get_spinor_from_field(in, spinprev, prev, timepos);
      
      get_su3matrix(&u,gaugefield,spacepos,timepos,3);
      get_su3matrix(&udagger,gaugefield,prev,timepos,3);
      adjoin_su3matrix(&udagger);

      if(coord[3] == NSPACE-1) spinor_apply_bc(spinnext, theta);
      else if(coord[3] == 0) spinor_apply_bc(spinprev, theta);
      
      multiply_spinor_gamma3(spinnext,tmp);
      real_multiply_spinor(tmp,-hmc_one_f);
      spinors_accumulate(spinnext,tmp);
      su3matrix_times_spinor(&u,spinnext,tmp);
      spinors_accumulate(spinout,tmp);

      multiply_spinor_gamma3(spinprev,tmp);
      spinors_accumulate(spinprev,tmp);
      su3matrix_times_spinor(&u,spinprev,tmp);
      spinors_accumulate(spinout,tmp);

      
      
      
      
      
      
      
      
      
      
      put_spinor_to_field(spinout,out,spacepos,timepos);
    }
  }

  return HMC_SUCCESS;
}

hmc_error dslash_eoprec(hmc_eoprec_spinor_field* out, hmc_eoprec_spinor_field* in, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im, int evenodd){
  
  for(int n=0; n<VOL4D/2; n++) {


    hmc_spinor spinout[SPINORSIZE];
    hmc_spinor spinnext[SPINORSIZE];
    hmc_spinor spinprev[SPINORSIZE];
    hmc_spinor tmp[SPINORSIZE];    
    hmc_su3matrix u;
    hmc_su3matrix udagger;
    int ns;
    int nt;
    int next;
    int prev;
    int neo_next;
    int neo_prev;
    
    set_local_zero_spinor(spinout);    

    //     int ns = get_nspace_from_eoprecindex(n,evenodd);
    //     int nt = get_ntime_from_eoprecindex(n,evenodd);
    if(evenodd == ODD) get_odd_site(n, &ns, &nt);
    else get_even_site(n, &ns, &nt);
       
    //like in host_geometry
    int coord[NDIM];
    coord[0]=0;
    for(int j=1;j<NDIM;j++) coord[j] = get_spacecoord(ns,j);    
    
    // spinout = U_0*(r-gamma_0)*spinnext + U^dagger_0(x-hat0) * (r+gamma_0)*spinprev
    next = (nt+1)%NTIME;
    prev = (nt-1+NTIME)%NTIME;
    
    neo_next = get_n_eoprec(next,ns);
    neo_prev = get_n_eoprec(prev,ns);

    get_spinor_from_eoprec_field(in,spinnext,neo_next);
    get_spinor_from_eoprec_field(in,spinprev,neo_prev);

    if(next == 0) spinor_apply_bc(spinnext, theta);
    else if(prev == NTIME-1) spinor_apply_bc(spinprev, theta);    
    
    get_su3matrix(&u,gaugefield,ns,nt,0);
    get_su3matrix(&udagger,gaugefield,ns,prev,0);
    adjoin_su3matrix(&udagger);
    
    //update links with chemical potential, this shall be put into compiler option lateron
    gaugefield_apply_chem_pot(&u, &udagger, chem_pot_re, chem_pot_im);
      
    multiply_spinor_gamma0(spinnext,tmp);
      real_multiply_spinor(tmp,-hmc_one_f);
      spinors_accumulate(spinnext,tmp);
      su3matrix_times_spinor(&u,spinnext,tmp);
      spinors_accumulate(spinout,tmp);

      multiply_spinor_gamma0(spinprev,tmp);
      spinors_accumulate(spinprev,tmp);
      su3matrix_times_spinor(&udagger,spinprev,tmp);
      spinors_accumulate(spinout,tmp);
     
    /*
    spinprojectproduct_gamma0(&u,spinnext,-hmc_one_f);
    spinors_accumulate(spinout,spinnext);

    spinprojectproduct_gamma0(&udagger,spinprev,hmc_one_f);
    spinors_accumulate(spinout,spinprev);
    */
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

    if(coord[1] == NSPACE-1) spinor_apply_bc(spinnext, theta);
    else if(coord[1] == 0) spinor_apply_bc(spinprev, theta);
     
  multiply_spinor_gamma1(spinnext,tmp);
      real_multiply_spinor(tmp,-hmc_one_f);
      spinors_accumulate(spinnext,tmp);
      su3matrix_times_spinor(&u,spinnext,tmp);
      spinors_accumulate(spinout,tmp);

      multiply_spinor_gamma1(spinprev,tmp);
      spinors_accumulate(spinprev,tmp);
      su3matrix_times_spinor(&u,spinprev,tmp);
      spinors_accumulate(spinout,tmp);  
      
    /*
    spinprojectproduct_gamma1(&u,spinnext,-hmc_one_f);
    spinors_accumulate(spinout,spinnext);

    spinprojectproduct_gamma1(&udagger,spinprev,hmc_one_f);
    spinors_accumulate(spinout,spinprev);
    */
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
    
    if(coord[2] == NSPACE-1) spinor_apply_bc(spinnext, theta);
    else if(coord[2] == 0) spinor_apply_bc(spinprev, theta);
      
    multiply_spinor_gamma2(spinnext,tmp);
      real_multiply_spinor(tmp,-hmc_one_f);
      spinors_accumulate(spinnext,tmp);
      su3matrix_times_spinor(&u,spinnext,tmp);
      spinors_accumulate(spinout,tmp);

      multiply_spinor_gamma2(spinprev,tmp);
      spinors_accumulate(spinprev,tmp);
      su3matrix_times_spinor(&u,spinprev,tmp);
      spinors_accumulate(spinout,tmp);   
    
    /*
    spinprojectproduct_gamma2(&u,spinnext,-hmc_one_f);
    spinors_accumulate(spinout,spinnext);

    spinprojectproduct_gamma2(&udagger,spinprev,hmc_one_f);
    spinors_accumulate(spinout,spinprev);
    */
    
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
      
    if(coord[3] == NSPACE-1) spinor_apply_bc(spinnext, theta);
    else if(coord[3] == 0) spinor_apply_bc(spinprev, theta);
    
    multiply_spinor_gamma3(spinnext,tmp);
      real_multiply_spinor(tmp,-hmc_one_f);
      spinors_accumulate(spinnext,tmp);
      su3matrix_times_spinor(&u,spinnext,tmp);
      spinors_accumulate(spinout,tmp);

      multiply_spinor_gamma3(spinprev,tmp);
      spinors_accumulate(spinprev,tmp);
      su3matrix_times_spinor(&u,spinprev,tmp);
      spinors_accumulate(spinout,tmp);
    
    /*
    spinprojectproduct_gamma3(&u,spinnext,-hmc_one_f);
    spinors_accumulate(spinout,spinnext);

    spinprojectproduct_gamma3(&udagger,spinprev,hmc_one_f);
    spinors_accumulate(spinout,spinprev);
    */ 
    
    //CP: why this kappa??
    real_multiply_spinor(spinout,-kappa);
    put_spinor_to_eoprec_field(spinout,out,n);

  }

  return HMC_SUCCESS;
}


hmc_error Aee(hmc_eoprec_spinor_field* in, hmc_eoprec_spinor_field* out, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im){
  hmc_eoprec_spinor_field* spintmp1 = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field* spintmp2 = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field* spintmp3 = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field* spintmp4 = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_complex one = hmc_complex_one;
  set_zero_eoprec_spinor(spintmp3);
  
//CP: original 
//   dslash_eoprec(spintmp1,in,gaugefield,kappa, theta, chem_pot_re, chem_pot_re,ODD); // D_oe
//   M_inverse_sitediagonal(spintmp2,spintmp1,kappa,mu); // R_o^(-1)
//   dslash_eoprec(out,spintmp2,gaugefield,kappa, theta, chem_pot_re, chem_pot_re,EVEN); // D_eo
//   M_sitediagonal(spintmp1,in,kappa,mu); //R_e
//   copy_spinor_eoprec(out, spintmp3);
//   saxpy_eoprec(spintmp3, spintmp1, &one, out); 

//CP: testing
  copy_spinor_eoprec(in, spintmp1);
  copy_spinor_eoprec(in, spintmp2);
  copy_spinor_eoprec(in, spintmp3);
  copy_spinor_eoprec(in, spintmp4);
  dslash_eoprec(spintmp1,in,gaugefield,kappa, theta, chem_pot_re, chem_pot_im,ODD); // D_oe
  M_inverse_sitediagonal(spintmp2,spintmp1,kappa,mu); // R_o^(-1)
  dslash_eoprec(spintmp3,spintmp2,gaugefield,kappa, theta, chem_pot_re, chem_pot_im,EVEN); // D_eo
  M_sitediagonal(spintmp4,in,kappa,mu); //R_e
 
//   copy_spinor_eoprec(spintmp3, spintmp2);
  copy_spinor_eoprec(spintmp1, out);
  
  
  delete [] spintmp1;
  delete [] spintmp2;
  delete [] spintmp3;
  delete [] spintmp4;

  return HMC_SUCCESS;
}