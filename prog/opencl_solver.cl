
void M_diag(hmc_spinor_field* in, hmc_spinor_field* out, hmc_float kappa, hmc_float mu){
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
  return; 
}

void dslash(hmc_spinor_field* in, hmc_spinor_field* out,__global hmc_ocl_gaugefield* gaugefield, hmc_float theta){

  hmc_spinor spinout[SPINORSIZE];

  for(int spacepos=0; spacepos<VOLSPACE; spacepos++) {
    for(int timepos=0; timepos<NTIME; timepos++) {

      hmc_spinor spinnext[SPINORSIZE];
      hmc_spinor spinprev[SPINORSIZE];
      hmc_spinor tmp[SPINORSIZE];
      hmc_ocl_su3matrix u;
      hmc_ocl_su3matrix udagger;
      int next;
      int prev;
      hmc_float theta = 0.;

      //like in host_geometry
      int coord[NDIM];
      coord[0]=0;
      for(int j=1;j<NDIM;j++) coord[j] = get_spacecoord(spacepos,j);
      
      set_local_zero_spinor(spinout);    
   
      // spinout = U_0*(r-gamma_0)*spinnext + U^dagger_0(x-hat0) * (r+gamma_0)*spinprev
      next = (timepos+1)%NTIME;
      prev = (timepos-1+NTIME)%NTIME;

      get_spinor_from_field(in, spinnext, spacepos, next);
      get_spinor_from_field(in, spinprev, spacepos, prev);

      if(next == 1) spinor_apply_bc(spinnext, theta);
      else if(prev == NTIME) spinor_apply_bc(spinprev, theta);
      
      get_su3matrix(&u,gaugefield,spacepos,timepos,0);
      get_su3matrix(&udagger,gaugefield,spacepos,prev,0);
      adjoin_su3matrix(&udagger);

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

      if(coord[1] == NSPACE) spinor_apply_bc(spinnext, theta);
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
      
      if(coord[2] == NSPACE) spinor_apply_bc(spinnext, theta);
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

      if(coord[3] == NSPACE) spinor_apply_bc(spinnext, theta);
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

  return;
}



void M(__private hmc_spinor_field* in, hmc_spinor_field* out,__global hmc_ocl_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta){
  
  M_diag(in, out, kappa, mu);    
  hmc_spinor_field tmp[SPINORFIELDSIZE];
  dslash(in,tmp,gaugefield, theta);

  hmc_complex kappa_cmplx = {kappa, 0.};
  saxpy(tmp, out, &kappa_cmplx, out);

  return;
}


void bicgstab(__global hmc_spinor_field* inout, hmc_eoprec_spinor_field* source, __global hmc_ocl_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta, int cgmax){

  //BiCGStab according to hep-lat/9404013

  hmc_spinor_field rn[SPINORFIELDSIZE];
  hmc_spinor_field rhat [SPINORFIELDSIZE];
  hmc_spinor_field v [SPINORFIELDSIZE];
  hmc_spinor_field p [SPINORFIELDSIZE];
  hmc_spinor_field s [SPINORFIELDSIZE];
  hmc_spinor_field t [SPINORFIELDSIZE];
  hmc_spinor_field inout_tmp [SPINORFIELDSIZE];
  hmc_complex rho;
  hmc_complex rho_next;
  hmc_complex alpha;
  hmc_complex omega;
  hmc_complex beta;

  hmc_complex tmp1;
  hmc_complex tmp2;
  hmc_complex one = hmc_complex_one;
  hmc_complex minusone = hmc_complex_minusone;
  
  //CP: copy spinor from global to local mem, this musst be avoided in the end
  get_spinor(inout, inout_tmp);

  for(int iter=0; iter<cgmax; iter++){

    if(iter%iter_refresh==0) {
      //fresh start
      M(inout_tmp,rn,gaugefield,kappa,mu, theta);
      saxpy(rn, source, &one, rn);
      copy_spinor(rn, rhat);

      alpha = hmc_complex_one;
      omega = hmc_complex_one;
      rho = hmc_complex_one;
      set_zero_spinor(v);
      set_zero_spinor(p);
    }

    rho_next = scalar_product(rhat,rn);
    tmp1 = complexdivide(&rho_next,&rho);
    rho = rho_next;
    tmp2 = complexdivide(&alpha,&omega);
    beta = complexmult(&tmp1,&tmp2);

    tmp1 = complexmult(&beta,&omega);
    tmp2 = complexmult(&minusone,&tmp1);
    saxsbypz(p, v, rn, &beta, &tmp2, p);

    M(p,v,gaugefield,kappa,mu, theta);

    tmp1 = scalar_product(rhat,v);
    alpha = complexdivide(&rho,&tmp1);

    saxpy(v, rn, &alpha, s);

    M(s,t,gaugefield,kappa,mu, theta);

    tmp1 = scalar_product(t,s);
    tmp2 = scalar_product(t,t);
    omega = complexdivide(&(tmp1),&(tmp2));

    saxpy(t, s, &omega, rn);

    saxsbypz(p, s, inout_tmp, &alpha, &omega, inout_tmp);

    hmc_float resid = global_squarenorm(rn);

    if(resid<epssquare) {
      hmc_spinor_field aux [SPINORFIELDSIZE];
      M(inout_tmp,aux,gaugefield,kappa,mu, theta);
      saxpy(aux, source, &one, aux);
      hmc_float trueresid = global_squarenorm(aux);

      if(trueresid<epssquare) return;
    }
  }
  
  //CP: copy spinor back to global mem, this musst be avoided in the end
  put_spinor(inout_tmp, inout);

  return;
}

void solver(__global hmc_spinor_field* in,__global hmc_spinor_field* out,__private hmc_spinor_field* b,__global hmc_ocl_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta, int cgmax){

  convert_to_kappa_format(in, kappa);
  bicgstab(in, b, gaugefield, kappa, mu, theta, cgmax);
  convert_from_kappa_format(in, out, kappa);
  
  return;
}


