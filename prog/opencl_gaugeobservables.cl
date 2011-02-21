
__kernel void plaquette(__global hmc_ocl_gaugefield * field,__global hmc_float * plaq_out, __global hmc_float* tplaq, __global hmc_float* splaq){
  int t, pos;
  int id = get_global_id(0);

  // this is an ugly workaround 'cuz there is no atomic_add for floats
  // FIXME replace by proper parallel reduction
  if( id > 0 )
	return;
  (*plaq_out) = 0.0f;
  (*splaq) = 0.0f;
  (*tplaq) = 0.0f;
  for( id = 0; id < VOL4D/2; ++id ) {

  hmc_float plaq=0;
  hmc_float splaq_tmp=0;
  hmc_float tplaq_tmp=0;
  hmc_float tmpfloat = 0;

  hmc_ocl_su3matrix tmp[SU3SIZE]; 
  hmc_ocl_su3matrix prod[SU3SIZE];
  
  //this is copied from the host-code (see comments there)
  //calc even plaquette
  get_even_site(id, &pos, &t);
  for(int mu=0; mu<NDIM; mu++) {
    for(int nu=0;nu<mu; nu++) {
      get_su3matrix(prod,field,pos,t,mu);
      if(mu==0) {
	int newt = (t+1)%NTIME;
	get_su3matrix(tmp,field,pos,newt,nu);
	} else {
	    get_su3matrix(tmp,field,get_neighbor(pos,mu),t,nu);
	}
	accumulate_su3matrix_prod(prod,tmp);
	if(nu==0) {
	  int newt = (t+1)%NTIME;
	  get_su3matrix(tmp,field,pos,newt,mu);
	} else {
	  get_su3matrix(tmp,field,get_neighbor(pos,nu),t,mu);
	}
	adjoin_su3matrix(tmp);
	accumulate_su3matrix_prod(prod,tmp);
	get_su3matrix(tmp,field,pos,t,nu);
	adjoin_su3matrix(tmp);
	accumulate_su3matrix_prod(prod,tmp);
	tmpfloat = trace_su3matrix(prod).re;
	plaq += tmpfloat;
	if(mu==0 || nu==0) {
	  tplaq_tmp+=tmpfloat;
	} else {
	  splaq_tmp+=tmpfloat;
	}
  }}
    
  //calc odd plaquette
  get_odd_site(id, &pos, &t);
  for(int mu=0; mu<NDIM; mu++) {
    for(int nu=0;nu<mu; nu++) {
      get_su3matrix(prod,field,pos,t,mu);
      if(mu==0) {
	int newt = (t+1)%NTIME;
	get_su3matrix(tmp,field,pos,newt,nu);
	} else {
	    get_su3matrix(tmp,field,get_neighbor(pos,mu),t,nu);
	}
	accumulate_su3matrix_prod(prod,tmp);
	if(nu==0) {
	  int newt = (t+1)%NTIME;
	  get_su3matrix(tmp,field,pos,newt,mu);
	} else {
	  get_su3matrix(tmp,field,get_neighbor(pos,nu),t,mu);
	}
	adjoin_su3matrix(tmp);
	accumulate_su3matrix_prod(prod,tmp);
	get_su3matrix(tmp,field,pos,t,nu);
	adjoin_su3matrix(tmp);
	accumulate_su3matrix_prod(prod,tmp);
	tmpfloat = trace_su3matrix(prod).re;
	plaq += tmpfloat;
	if(mu==0 || nu==0) {
	  tplaq_tmp+=tmpfloat;
	} else {
	  splaq_tmp+=tmpfloat;
	}
  }}
  
  (*plaq_out) += plaq;
  (*splaq) += splaq_tmp;
  (*tplaq) += tplaq_tmp;

  }
}


__kernel void polyakov(__global hmc_ocl_gaugefield * field, __global hmc_complex * out){
  int const tdir = 0;
  int pos, t;
  int id = get_global_id(0);
   
  // this is an ugly workaround 'cuz there is no atomic_add for floats
  // FIXME replace by proper parallel reduction
  if( id > 0 )
	return;
  (*out).re = 0.0f;
  (*out).im = 0.0f;
  for( id = 0; id < VOL4D/2; ++id ) {

  hmc_ocl_su3matrix prod[SU3SIZE];
  hmc_ocl_su3matrix tmp[SU3SIZE];
  unit_su3matrix(prod);
  hmc_complex tmpcomplex;
  
  get_even_site(id, &pos, &t);
  if(t==0){
    for(int t=0; t<NTIME; t++) {
      get_su3matrix(tmp,field,pos,t,tdir);
      accumulate_su3matrix_prod(prod,tmp);
    }
    tmpcomplex = trace_su3matrix(prod);
    (*out).re += tmpcomplex.re;
    (*out).im += tmpcomplex.im;
  }
  else if(t==1){
    get_odd_site(id, &pos, &t);
    for(int t=0; t<NTIME; t++) {
      get_su3matrix(tmp,field,pos,t,tdir);
      accumulate_su3matrix_prod(prod,tmp);
    }
    tmpcomplex = trace_su3matrix(prod);
    (*out).re += tmpcomplex.re;
    (*out).im += tmpcomplex.im;
  }
  else continue; //return;

  }
}
