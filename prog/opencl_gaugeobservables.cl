//opencl_gaugeobservables.cl

__kernel void plaquette(__global hmc_ocl_gaugefield * field,__global hmc_float * plaq_out, __global hmc_float* tplaq, __global hmc_float* splaq){

  int t, pos, id, id_tmp, size;
  id_tmp = get_global_id(0);
  size = get_global_size(0);

  // this is an ugly workaround 'cuz there is no atomic_add for floats
  // FIXME replace by proper parallel reduction
  if( id > 0 )
	return;
  (*plaq_out) = 0.0f;
  (*splaq) = 0.0f;
  (*tplaq) = 0.0f;
  for( id = 0; id < get_global_size(0); ++id ) {

  hmc_float plaq=0;
  hmc_float splaq_tmp=0;
  hmc_float tplaq_tmp=0;
  hmc_float tmpfloat = 0;

  hmc_ocl_su3matrix tmp[SU3SIZE]; 
  hmc_ocl_su3matrix prod[SU3SIZE];

  for(id = id_tmp; id<VOLSPACE*NTIME/2; id+=size){
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
  }
  
// TODO use reduction
//  (plaq_out)[id_tmp] += plaq;
//  (splaq)[id_tmp] += splaq_tmp;
//  (tplaq)[id_tmp] += tplaq_tmp;
//  
//  //wait for all threads to end calculations, does this work in a kernel???
//  //cl_finish(queue);
//  
//  //perform reduction
//  int cut1;
//  int cut2 = size;
//  if(size > 128){
//    for(cut1 = 128; cut1>0; cut1/=2){
//      for(int i = id_tmp+cut1; i < cut2; i+=cut1){
//	(plaq_out)[id_tmp] += (plaq_out)[i];
//	(splaq)[id_tmp] += (splaq)[i];
//	(tplaq)[id_tmp] += (tplaq)[i];
//      }
//      cut2 = cut1;
//    }
//  }
//  else if(id_tmp == 0) {
//    for(int i = id_tmp; i < size; i++){
//      (plaq_out)[id_tmp] += (plaq_out)[i];
//      (splaq)[id_tmp] += (splaq)[i];
//      (tplaq)[id_tmp] += (tplaq)[i];
//    }
//  }
//    
//  return;
  (*plaq_out) += plaq;
  (*splaq) += splaq_tmp;
  (*tplaq) += tplaq_tmp;

  }
}

__kernel void polyakov(__global hmc_ocl_gaugefield * field, __global hmc_complex * out){
  
  int t, pos, id, size;
  id = get_global_id(0);
  size = get_global_size(0);
  int const tdir = 0;
   
  // this is an ugly workaround 'cuz there is no atomic_add for floats
  // FIXME replace by proper parallel reduction
  if( id > 0 )
	return;
  (*out).re = 0.0f;
  (*out).im = 0.0f;
  for( id = 0; id < VOLSPACE/2; ++id ) {

  hmc_ocl_su3matrix prod[SU3SIZE];
  hmc_ocl_su3matrix tmp[SU3SIZE];
  unit_su3matrix(prod);
  hmc_complex tmpcomplex;
  hmc_complex tmp_pol;
  tmp_pol.re = 0.;
  tmp_pol.im = 0.;
  
  //for(id = id_tmp; id<VOLSPACE/2; id+=size){
    
    //calc polyakov loop at even site
    get_even_site(id, &pos, &t);
    for(int t=0; t<NTIME; t++) {
      get_su3matrix(tmp,field,pos,t,tdir);
      accumulate_su3matrix_prod(prod,tmp);
    }
    tmpcomplex = trace_su3matrix(prod);
    (tmp_pol).re += tmpcomplex.re;
    (tmp_pol).im += tmpcomplex.im;
    
    //calc polyakov loop at odd site
    get_odd_site(id, &pos, &t);
    for(int t=0; t<NTIME; t++) {
      get_su3matrix(tmp,field,pos,t,tdir);
      accumulate_su3matrix_prod(prod,tmp);
    }
    tmpcomplex = trace_su3matrix(prod);
    (tmp_pol).re += tmpcomplex.re;
    (tmp_pol).im += tmpcomplex.im;
  }
  ((out)[id_tmp]).re += tmp_pol.re;
  ((out)[id_tmp]).im += tmp_pol.im;
  
// TODO use reduction
//  //wait for all threads to end calculations, does this work in a kernel???
//  //cl_finish(queue);
//  
//  //perform reduction
//  int cut1;
//  int cut2 = size;
//  if(size > 128){
//    for(cut1 = 128; cut1>0; cut1/=2){
//      for(int i = id_tmp+cut1; i < cut2; i+=cut1){
//	((out)[id_tmp]).re +=  ((out)[i]).re;
//	((out)[id_tmp]).im +=  ((out)[i]).im;
//      }
//      cut2 = cut1;
//    }
//  }
//  else if(id_tmp == 0) {
//    for(int i = id_tmp; i < size; i++){
//     ((out)[id_tmp]).re +=  ((out)[i]).re;
//     ((out)[id_tmp]).im +=  ((out)[i]).im;
//    }
//  }
//  
//  return;
    else continue; //return;

  }
}
