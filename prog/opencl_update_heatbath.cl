

void calc_staple(__global hmc_ocl_gaugefield* field,__private hmc_ocl_staplematrix* dest, const int pos, const int t, const int mu_in){
  hmc_ocl_su3matrix prod[SU3SIZE];
  hmc_ocl_su3matrix prod2[SU3SIZE];
  hmc_ocl_su3matrix tmp[SU3SIZE];
  hmc_ocl_staplematrix dummy[STAPLEMATRIXSIZE];
  int nu, newpos, newt;
  
  zero_staplematrix(dummy);
  
  //iterate through the three directions other than mu
  for(int i = 1; i<NDIM; i++){
    
    zero_su3matrix(prod);
    nu = (mu_in + i)%NDIM;
    //first staple
    //u_nu(x+mu)
    if(mu_in==0) {
      newt = (t+1)%NTIME;
      get_su3matrix(tmp,field,pos,newt,nu);
    } else {
      get_su3matrix(tmp,field,get_neighbor(pos,mu_in),t,nu);
    }
//     get_su3matrix(tmp,field,1,0,0);
    
    copy_su3matrix(prod, tmp);
    //adjoint(u_mu(x+nu))
    if(nu==0) {
      newt = (t+1)%NTIME;
      get_su3matrix(tmp,field,pos,newt,mu_in);
    } else {
      get_su3matrix(tmp,field,get_neighbor(pos,nu),t,mu_in);
    }
    adjoin_su3matrix(tmp);
    accumulate_su3matrix_prod(prod,tmp);
    //adjoint(u_nu(x))
    get_su3matrix(tmp,field,pos,t,nu);
    adjoin_su3matrix(tmp);
    accumulate_su3matrix_prod(prod,tmp);  
    //second staple
    //adjoint (u_nu(x+mu-nu))
    //newpos is "pos-nu" (spatial)
    newpos = get_lower_neighbor(pos, nu);
    if(mu_in==0) {
      newt = (t+1)%NTIME;
      get_su3matrix(tmp,field,newpos,newt,nu);
    } else if (nu == 0){
        newt = (t-1+NTIME)%NTIME;
        get_su3matrix(tmp,field,get_neighbor(pos,mu_in),newt,nu);
    }
    else{
	get_su3matrix(tmp,field,get_neighbor(newpos,mu_in),t,nu);
    } 
    adjoin_su3matrix(tmp);
    copy_su3matrix(prod2, tmp);
    //adjoint(u_mu(x-nu))
    if(mu_in==0) {
      get_su3matrix(tmp,field,newpos,t,mu_in);
    } 
    else if (nu == 0){
        newt = (t-1+NTIME)%NTIME;
        get_su3matrix(tmp,field,pos,newt,mu_in);
    }
    else{
	get_su3matrix(tmp,field,newpos,t,mu_in);
    }
    adjoin_su3matrix(tmp);
    accumulate_su3matrix_prod(prod2,tmp);
    //adjoint(u_nu(x-nu))
    if(mu_in==0) {
      get_su3matrix(tmp,field,newpos,t,nu);
    } 
    else if (nu == 0){
        newt = (t-1+NTIME)%NTIME;
        get_su3matrix(tmp,field,pos,newt,nu);
    }
    else{
	get_su3matrix(tmp,field,newpos,t,nu);
    }
    accumulate_su3matrix_prod(prod2,tmp); 
   
    accumulate_su3matrices_add(dummy, prod);
    accumulate_su3matrices_add(dummy, prod2);
  }
  copy_staplematrix(dest, dummy);
}

__kernel void heatbath_even(__global hmc_ocl_gaugefield* gaugefield, const hmc_float beta, const int mu, __global hmc_ocl_ran * rnd){
  //it is assumed that there are VOL4D/2 threads
  int t, pos;
  int id = get_global_id(0);
  get_even_site(id, &pos, &t);
  
  hmc_ocl_su3matrix U[SU3SIZE];
  hmc_ocl_staplematrix W[STAPLEMATRIXSIZE];
  hmc_ocl_staplematrix staplematrix[STAPLEMATRIXSIZE];
  
  hmc_complex w [su2_entries];
  hmc_float w_pauli[su2_entries], k, r_pauli[su2_entries];
  int order[3]; 
  random_1_2_3(order, &rnd[id]);
  get_su3matrix(U, gaugefield, pos, t, mu);
  
  hmc_complex det = det_su3matrix(U);
  hmc_complex detadj = complexconj(&det);
  hmc_complex detsqnorm = complexmult(&det, &detadj);
  if( (detsqnorm.re - hmc_one_f) <= projectioneps)
    project_su3(U);

  calc_staple(gaugefield, staplematrix, pos, t, mu);

  for(int i=0; i<3; i++)
  {
    multiply_staplematrix(W, U, staplematrix); 
    reduction(w, W, order[i]);

    w_pauli[0] = 0.5*(w[0].re + w[3].re);
    w_pauli[1] = 0.5*(w[1].im + w[2].im);
    w_pauli[2] = 0.5*(w[1].re - w[2].re);
    w_pauli[3] = 0.5*(w[0].im - w[3].im);
    k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );

    hmc_float beta_neu =  2.*beta / NC*k;
    SU2Update(r_pauli, beta_neu, &rnd[id]);
    
    w[0].re = (r_pauli[0]*w_pauli[0] + r_pauli[1]*w_pauli[1] + r_pauli[2]*w_pauli[2] + r_pauli[3]*w_pauli[3] )/k;
    w[0].im = (w_pauli[0]*r_pauli[3] - w_pauli[3]*r_pauli[0] + r_pauli[1]*w_pauli[2] - r_pauli[2]*w_pauli[1] )/k;
    w[1].re = (w_pauli[0]*r_pauli[2] - w_pauli[2]*r_pauli[0] + r_pauli[3]*w_pauli[1] - r_pauli[1]*w_pauli[3] )/k;
    w[1].im = (w_pauli[0]*r_pauli[1] - w_pauli[1]*r_pauli[0] + r_pauli[2]*w_pauli[3] - r_pauli[3]*w_pauli[2] )/k;
    w[2].re = -(w_pauli[0]*r_pauli[2] - w_pauli[2]*r_pauli[0] + r_pauli[3]*w_pauli[1] - r_pauli[1]*w_pauli[3] )/k;
    w[2].im = (w_pauli[0]*r_pauli[1] - w_pauli[1]*r_pauli[0] + r_pauli[2]*w_pauli[3] - r_pauli[3]*w_pauli[2] )/k;
    w[3].re = (r_pauli[0]*w_pauli[0] + r_pauli[1]*w_pauli[1] + r_pauli[2]*w_pauli[2] + r_pauli[3]*w_pauli[3] )/k;
    w[3].im = -(w_pauli[0]*r_pauli[3] - w_pauli[3]*r_pauli[0] + r_pauli[1]*w_pauli[2] - r_pauli[2]*w_pauli[1] )/k;
    
    hmc_ocl_su3matrix extW[SU3SIZE]; 
    extend (extW, order[i], w); 
    //perhaps one can define another acc-func here and save one copying step
    accumulate_su3matrix_prod(extW, U);
    copy_su3matrix(U, extW);
  }
  put_su3matrix(gaugefield, U, pos, t, mu);

  return;
}

__kernel void heatbath_odd(__global hmc_ocl_gaugefield* gaugefield, const hmc_float beta, const int mu, __global hmc_ocl_ran * rnd){

  int t, pos;
  int id = get_global_id(0);
  get_odd_site(id, &pos, &t);
  
  hmc_ocl_su3matrix U[SU3SIZE];
  hmc_ocl_staplematrix W[STAPLEMATRIXSIZE];
  hmc_ocl_staplematrix staplematrix[STAPLEMATRIXSIZE];
  
  hmc_complex w [su2_entries];
  hmc_float w_pauli[su2_entries], k, r_pauli[su2_entries];
  int order[3]; 
  random_1_2_3(order, &rnd[id]);
  //old U
  get_su3matrix(U, gaugefield, pos, t, mu);
     
  hmc_complex det = det_su3matrix(U);
  hmc_complex detadj = complexconj(&det);
  hmc_complex detsqnorm = complexmult(&det, &detadj);
  if( (detsqnorm.re - hmc_one_f) <= projectioneps)
    project_su3(U);

  calc_staple(gaugefield, staplematrix, pos, t, mu);

  for(int i=0; i<3; i++)
  {
    multiply_staplematrix(W, U, staplematrix); 
    reduction(w, W, order[i]);

    w_pauli[0] = 0.5*(w[0].re + w[3].re);
    w_pauli[1] = 0.5*(w[1].im + w[2].im);
    w_pauli[2] = 0.5*(w[1].re - w[2].re);
    w_pauli[3] = 0.5*(w[0].im - w[3].im);
    k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );

    hmc_float beta_neu =  2.*beta / NC*k;
    SU2Update(r_pauli, beta_neu, &rnd[id]);
    
    w[0].re = (r_pauli[0]*w_pauli[0] + r_pauli[1]*w_pauli[1] + r_pauli[2]*w_pauli[2] + r_pauli[3]*w_pauli[3] )/k;
    w[0].im = (w_pauli[0]*r_pauli[3] - w_pauli[3]*r_pauli[0] + r_pauli[1]*w_pauli[2] - r_pauli[2]*w_pauli[1] )/k;
    w[1].re = (w_pauli[0]*r_pauli[2] - w_pauli[2]*r_pauli[0] + r_pauli[3]*w_pauli[1] - r_pauli[1]*w_pauli[3] )/k;
    w[1].im = (w_pauli[0]*r_pauli[1] - w_pauli[1]*r_pauli[0] + r_pauli[2]*w_pauli[3] - r_pauli[3]*w_pauli[2] )/k;
    w[2].re = -(w_pauli[0]*r_pauli[2] - w_pauli[2]*r_pauli[0] + r_pauli[3]*w_pauli[1] - r_pauli[1]*w_pauli[3] )/k;
    w[2].im = (w_pauli[0]*r_pauli[1] - w_pauli[1]*r_pauli[0] + r_pauli[2]*w_pauli[3] - r_pauli[3]*w_pauli[2] )/k;
    w[3].re = (r_pauli[0]*w_pauli[0] + r_pauli[1]*w_pauli[1] + r_pauli[2]*w_pauli[2] + r_pauli[3]*w_pauli[3] )/k;
    w[3].im = -(w_pauli[0]*r_pauli[3] - w_pauli[3]*r_pauli[0] + r_pauli[1]*w_pauli[2] - r_pauli[2]*w_pauli[1] )/k;
    
    hmc_ocl_su3matrix extW[SU3SIZE]; 
    extend (extW, order[i], w); 
    //perhaps one can define another acc-func here and save one copying step
    accumulate_su3matrix_prod(extW, U);
    copy_su3matrix(U, extW);
  }
  put_su3matrix(gaugefield, U, pos, t, mu);

  return;
}

__kernel void overrelax_even(__global hmc_ocl_gaugefield* gaugefield, const hmc_float beta, const int mu, __global hmc_ocl_ran * rnd){
  //it is assumed that there are VOL4D/2 threads
  int t, pos;
  int id = get_global_id(0);
  get_even_site(id, &pos, &t);
  
  hmc_ocl_su3matrix U[SU3SIZE];
  hmc_ocl_staplematrix W[STAPLEMATRIXSIZE];
  hmc_ocl_staplematrix staplematrix[STAPLEMATRIXSIZE];
  
  hmc_complex w [su2_entries];
  hmc_float w_pauli[su2_entries], k;
  int order[3]; 
  random_1_2_3(order, &rnd[id]);
  get_su3matrix(U, gaugefield, pos, t, mu);

  hmc_complex det = det_su3matrix(U);
  hmc_complex detadj = complexconj(&det);
  hmc_complex detsqnorm = complexmult(&det, &detadj);
  if( (detsqnorm.re - hmc_one_f) <= projectioneps)
    project_su3(U);
  
  calc_staple(gaugefield, staplematrix, pos, t, mu);
  hmc_ocl_su3matrix tmp[SU3SIZE]; 
  hmc_ocl_su3matrix extW[SU3SIZE]; 
  for(int i=0; i<3; i++)
  {
    multiply_staplematrix(W, U, staplematrix); 
    reduction(w, W, order[i]);

    w_pauli[0] = 0.5*(w[0].re + w[3].re);
    w_pauli[1] = 0.5*(w[1].im + w[2].im);
    w_pauli[2] = 0.5*(w[1].re - w[2].re);
    w_pauli[3] = 0.5*(w[0].im - w[3].im);
    k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );

    w[0].re = (w_pauli[0]*w_pauli[0] - w_pauli[1]*w_pauli[1] - w_pauli[2]*w_pauli[2] - w_pauli[3]*w_pauli[3])/k/k;
    w[0].im = (-2.*w_pauli[0]*w_pauli[3])/k/k;
    w[1].re = (-2.*w_pauli[0]*w_pauli[2])/k/k;
    w[1].im = (-2.*w_pauli[0]*w_pauli[1])/k/k;
    w[2].re = (2.*w_pauli[0]*w_pauli[2])/k/k;
    w[2].im = (-2.*w_pauli[0]*w_pauli[1])/k/k;
    w[3].re = (w_pauli[0]*w_pauli[0] - w_pauli[1]*w_pauli[1] - w_pauli[2]*w_pauli[2] - w_pauli[3]*w_pauli[3])/k/k;
    w[3].im = (2.*w_pauli[0]*w_pauli[3])/k/k;
    
    extend (extW, order[i], w); 
    multiply_su3matrices(tmp, extW, U);
    copy_su3matrix(U, tmp);
  }
  put_su3matrix(gaugefield, U, pos, t, mu);
  return;
}

__kernel void overrelax_odd(__global hmc_ocl_gaugefield* gaugefield, const hmc_float beta, const int mu, __global hmc_ocl_ran * rnd){
  int t, pos;
  int id = get_global_id(0);
  get_odd_site(id, &pos, &t);
  
  hmc_ocl_su3matrix U[SU3SIZE];
  hmc_ocl_staplematrix W[STAPLEMATRIXSIZE];
  hmc_ocl_staplematrix staplematrix[STAPLEMATRIXSIZE];
  
  hmc_complex w [su2_entries];
  hmc_float w_pauli[su2_entries], k;
  int order[3]; 
  random_1_2_3(order, &rnd[id]);
  get_su3matrix(U, gaugefield, pos, t, mu);

  hmc_complex det = det_su3matrix(U);
  hmc_complex detadj = complexconj(&det);
  hmc_complex detsqnorm = complexmult(&det, &detadj);
  if( (detsqnorm.re - hmc_one_f) <= projectioneps)
    project_su3(U);
  
  calc_staple(gaugefield, staplematrix, pos, t, mu);
  hmc_ocl_su3matrix tmp[SU3SIZE]; 
  hmc_ocl_su3matrix extW[SU3SIZE]; 
  for(int i=0; i<3; i++)
  {
    multiply_staplematrix(W, U, staplematrix); 
    reduction(w, W, order[i]);

    w_pauli[0] = 0.5*(w[0].re + w[3].re);
    w_pauli[1] = 0.5*(w[1].im + w[2].im);
    w_pauli[2] = 0.5*(w[1].re - w[2].re);
    w_pauli[3] = 0.5*(w[0].im - w[3].im);
    k = sqrt(  w_pauli[0]*w_pauli[0] +  w_pauli[1]*w_pauli[1] + w_pauli[2]*w_pauli[2] + w_pauli[3]*w_pauli[3]  );

    w[0].re = (w_pauli[0]*w_pauli[0] - w_pauli[1]*w_pauli[1] - w_pauli[2]*w_pauli[2] - w_pauli[3]*w_pauli[3])/k/k;
    w[0].im = (-2.*w_pauli[0]*w_pauli[3])/k/k;
    w[1].re = (-2.*w_pauli[0]*w_pauli[2])/k/k;
    w[1].im = (-2.*w_pauli[0]*w_pauli[1])/k/k;
    w[2].re = (2.*w_pauli[0]*w_pauli[2])/k/k;
    w[2].im = (-2.*w_pauli[0]*w_pauli[1])/k/k;
    w[3].re = (w_pauli[0]*w_pauli[0] - w_pauli[1]*w_pauli[1] - w_pauli[2]*w_pauli[2] - w_pauli[3]*w_pauli[3])/k/k;
    w[3].im = (2.*w_pauli[0]*w_pauli[3])/k/k;
    
    extend (extW, order[i], w); 
    //perhaps one can define another acc-func here and save one copying step
    multiply_su3matrices(tmp, extW, U);
    copy_su3matrix(U, tmp);
  }
  put_su3matrix(gaugefield, U, pos, t, mu);
  return;
}