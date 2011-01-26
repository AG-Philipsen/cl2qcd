
//CP: define new functions to speed up the alg



get_adjsu3matrix
//out = adj(in)
//to be done
void copy_adjsu3matrix(__private hmc_ocl_su3matrix * out, __private hmc_ocl_su3matrix * in){
  for(int n = 0; n<SU3SIZE; n++) {
    (out[n]).re = (in[n]).re;
    (out[n]).im = (in[n]).im;
  }
  return;
}

//U = U*adj(V)
//to be done
void accumulate_su3matrix_adjprod(__private hmc_ocl_su3matrix *acc,__private hmc_ocl_su3matrix *mult){
  hmc_ocl_su3matrix tmp[SU3SIZE];
  multiply_adjsu3matrices(tmp,acc,mult);
  copy_su3matrix(acc,tmp);
  return;
}

void calc_staple(__global hmc_ocl_gaugefield* field,__private hmc_ocl_staplematrix* dest, const int pos, const int t, const int mu_in){
  
  //assume that dest is set to zero!!!!!
  
  hmc_ocl_su3matrix prod[SU3SIZE];
  hmc_ocl_su3matrix prod2[SU3SIZE];
  hmc_ocl_su3matrix tmp[SU3SIZE];
  hmc_ocl_su3matrix tmp2[SU3SIZE];
  int nu, newpos, newt;
  
  zero_staplematrix(dummy);
  
  //iterate through the three directions other than mu
  for(int i = 1; i<NDIM; i++){
    zero_su3matrix(prod);
    nu = (mu_in + i)%NDIM;
    newpos = get_lower_neighbor(pos, nu);
    if(mu_in==0) {
      newt = (t+1)%NTIME;
      get_su3matrix(prod,field,pos,newt,nu);
      get_adjsu3matrix(prod2,field,newpos,newt,nu);
      get_su3matrix(tmp,field,get_neighbor(pos,nu),t,mu_in);
      get_su3matrix(tmp2,field,newpos,t,mu_in);
      accumulate_su3matrix_adjprod(prod,tmp);
      accumulate_su3matrix_adjprod(prod2,tmp2); 
      get_su3matrix(tmp,field,pos,t,nu);
      get_su3matrix(tmp2,field,newpos,t,nu);
      accumulate_su3matrix_adjprod(prod,tmp);  
      accumulate_su3matrix_adjprod(prod2,tmp2); 
    } 
    else if(nu == 0){
      get_su3matrix(prod,field,get_neighbor(pos,mu_in),t,nu);
      newt = (t-1+NTIME)%NTIME;
      get_adjsu3matrix(prod2,field,get_neighbor(pos,mu_in),newt,nu);
      newt = (t+1)%NTIME;
      get_su3matrix(tmp,field,pos,newt,mu_in);
      newt = (t-1+NTIME)%NTIME;
      get_su3matrix(tmp2,field,pos,newt,mu_in);
      accumulate_su3matrix_adjprod(prod,tmp);
      accumulate_su3matrix_adjprod(prod2,tmp2); 
      get_su3matrix(tmp,field,pos,t,nu);
      newt = (t-1+NTIME)%NTIME;
      get_su3matrix(tmp2,field,pos,newt,nu);
      accumulate_su3matrix_adjprod(prod,tmp);  
      accumulate_su3matrix_adjprod(prod2,tmp2); 
    }
    else{
      get_su3matrix(prod,field,get_neighbor(pos,mu_in),t,nu);
      get_adjsu3matrix(prod2,field,get_neighbor(newpos,mu_in),t,nu);
      get_su3matrix(tmp,field,get_neighbor(pos,nu),t,mu_in);
      get_su3matrix(tmp2,field,newpos,t,mu_in);
      accumulate_su3matrix_adjprod(prod,tmp);
      accumulate_su3matrix_adjprod(prod2,tmp2); 
      get_su3matrix(tmp,field,pos,t,nu);
      get_su3matrix(tmp2,field,newpos,t,nu);
      accumulate_su3matrix_adjprod(prod,tmp);  
      accumulate_su3matrix_adjprod(prod2,tmp2); 
    }
    accumulate_su3matrices_add(dest, prod);
    accumulate_su3matrices_add(dest, prod2);
  }
}