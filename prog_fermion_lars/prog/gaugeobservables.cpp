#include "gaugeobservables.h"
#include <cstdio>

void print_gaugeobservables(hmc_gaugefield* field){
  hmc_float tplaq=0;
  hmc_float splaq=0;
  hmc_float plaq = plaquette(field,&tplaq,&splaq);
  hmc_complex pol = polyakov(field);
  printf("%f\t%f\t%f\t%f\t%f\n",plaq,tplaq,splaq,pol.re,pol.im);
  return;
}

hmc_float plaquette(hmc_gaugefield * field){
  hmc_float tdummy;
  hmc_float sdummy;
  return plaquette(field,&tdummy,&sdummy);
}

hmc_float plaquette(hmc_gaugefield * field, hmc_float* tplaq, hmc_float* splaq){
  hmc_float plaq=0;
  *tplaq=0;
  *splaq=0;
  for(int t=0;t<NTIME;t++) {
    for(int n=0;n<VOLSPACE;n++) {
      for(int mu=0; mu<NDIM; mu++) {
	for(int nu=0;nu<mu; nu++) {
	  hmc_su3matrix tmp;
	  hmc_su3matrix prod;
	  //u_mu(x)
	  get_su3matrix(&prod,field,n,t,mu);
	  //u_nu(x+mu)
	  if(mu==0) {
	    int newt = (t+1)%NTIME; //(haha)
	    get_su3matrix(&tmp,field,n,newt,nu);
	  } else {
	    get_su3matrix(&tmp,field,get_neighbor(n,mu),t,nu);
	  }
	  accumulate_su3matrix_prod(&prod,&tmp);
	  //adjoint(u_mu(x+nu))
	  if(nu==0) {
	    int newt = (t+1)%NTIME;
	    get_su3matrix(&tmp,field,n,newt,mu);
	  } else {
	    get_su3matrix(&tmp,field,get_neighbor(n,nu),t,mu);
	  }
	  adjoin_su3matrix(&tmp);
	  accumulate_su3matrix_prod(&prod,&tmp);
	  //adjoint(u_nu(x))
	  get_su3matrix(&tmp,field,n,t,nu);
	  adjoin_su3matrix(&tmp);
	  accumulate_su3matrix_prod(&prod,&tmp);
	  hmc_float tmpfloat = trace_su3matrix(&prod).re;
	  plaq += tmpfloat;
	  if(mu==0 || nu==0) {
	    *tplaq+=tmpfloat;
	  } else {
	    *splaq+=tmpfloat;
	  }
	}
      }
    }
  }
  *tplaq /= static_cast<hmc_float>(VOL4D*NC*(NDIM-1));
  *splaq /= static_cast<hmc_float>(VOL4D*NC*(NDIM-1)*(NDIM-2))/2. ;
  return plaq*2.0/static_cast<hmc_float>(VOL4D*NDIM*(NDIM-1)*NC);
}


hmc_complex polyakov(hmc_gaugefield * field){
  int const tdir = 0;
  hmc_complex res;
  res.re=0;
  res.im=0;
  for(int n=0; n<VOLSPACE; n++) {
    hmc_su3matrix prod;
    unit_su3matrix(&prod);
    for(int t=0; t<NTIME; t++) {
      hmc_su3matrix tmp;
      get_su3matrix(&tmp,field,n,t,tdir);
      accumulate_su3matrix_prod(&prod,&tmp);
    }
    hmc_complex tmpcomplex = trace_su3matrix(&prod);
    complexaccumulate(&res,&tmpcomplex);
  }
  res.re /= static_cast<hmc_float>(NC*VOLSPACE);
  res.im /= static_cast<hmc_float>(NC*VOLSPACE);
  return res;
}


//!!not tested
hmc_complex spatial_polyakov(hmc_gaugefield* field, int dir){
  //assuming dir=1,2, or 3
  hmc_complex res;
  res.re=0;
  res.im=0;
  for(int x1=0; x1<NSPACE; x1++) {
    for(int x2=0; x2<NSPACE; x2++) {
      for(int t=0; t<NTIME; t++) {
	hmc_su3matrix prod;
	unit_su3matrix(&prod);
	for(int xpol=0; xpol<NSPACE; xpol++) {
	  hmc_su3matrix tmp;
	  int coord[NDIM];
	  coord[0]=t;
	  coord[dir]=xpol;
	  int next = (dir%(NDIM-1)) + 1;
	  coord[next]=x1;
	  int nnext = (next%(NDIM-1)) + 1;
	  coord[nnext]=x2;
	  int pos = get_nspace(coord);
	  get_su3matrix(&tmp,field,pos,t,dir);
	  accumulate_su3matrix_prod(&prod,&tmp);
	}
	hmc_complex tmpcomplex = trace_su3matrix(&prod);
	complexaccumulate(&res,&tmpcomplex);
      }
    }
  }
  res.re /= static_cast<hmc_float>(NC*NSPACE*NSPACE*NTIME);
  res.im /= static_cast<hmc_float>(NC*NSPACE*NSPACE*NTIME);
  return res;
}

//!! not tested
hmc_complex polyakov_x(hmc_gaugefield * field){
  int const tdir = 1; 
  int pos;
  hmc_complex res;
  res.re=0;
  res.im=0;
  for(int n=0; n<NSPACE*NSPACE; n++) {
    for(int t=0; t<NTIME; t++) {
    hmc_su3matrix prod;
    unit_su3matrix(&prod);
    for(int x=0; x<NSPACE; x++) {
      pos = x + NSPACE*n;
      hmc_su3matrix tmp;
      get_su3matrix(&tmp,field,pos,t,tdir);
      accumulate_su3matrix_prod(&prod,&tmp);
    }
    hmc_complex tmpcomplex = trace_su3matrix(&prod);
    complexaccumulate(&res,&tmpcomplex);
  }}
  res.re /= static_cast<hmc_float>(NC*NSPACE*NSPACE*NTIME);
  res.im /= static_cast<hmc_float>(NC*NSPACE*NSPACE*NTIME);
  return res;
}

//!! not tested
hmc_complex polyakov_y(hmc_gaugefield * field){
  int const tdir = 2; 
  int pos;
  hmc_complex res;
  res.re=0;
  res.im=0;
  for(int x=0; x<NSPACE; x++) {
    for(int z=0; z<NSPACE; z++) {
    for(int t=0; t<NTIME; t++) {
    hmc_su3matrix prod;
    unit_su3matrix(&prod);
    for(int y=0; y<NSPACE; y++) {
      pos = x + y*NSPACE + z*NSPACE*NSPACE;
      hmc_su3matrix tmp;
      get_su3matrix(&tmp,field,pos,t,tdir);
      accumulate_su3matrix_prod(&prod,&tmp);
    }
    hmc_complex tmpcomplex = trace_su3matrix(&prod);
    complexaccumulate(&res,&tmpcomplex);
  }}}
  res.re /= static_cast<hmc_float>(NC*NSPACE*NSPACE*NTIME);
  res.im /= static_cast<hmc_float>(NC*NSPACE*NSPACE*NTIME);
  return res;
}

//!! not tested
hmc_complex polyakov_z(hmc_gaugefield * field){
  int const tdir = 3; 
  int pos;
  hmc_complex res;
  res.re=0;
  res.im=0;
  for(int n=0; n<NSPACE*NSPACE; n++) {
    for(int t=0; t<NTIME; t++) {
    hmc_su3matrix prod;
    unit_su3matrix(&prod);
    for(int z=0; z<NSPACE; z++) {
      pos = n + z*NSPACE*NSPACE;
      hmc_su3matrix tmp;
      get_su3matrix(&tmp,field,pos,t,tdir);
      accumulate_su3matrix_prod(&prod,&tmp);
    }
    hmc_complex tmpcomplex = trace_su3matrix(&prod);
    complexaccumulate(&res,&tmpcomplex);
  }}
  res.re /= static_cast<hmc_float>(NC*NSPACE*NSPACE*NTIME);
  res.im /= static_cast<hmc_float>(NC*NSPACE*NSPACE*NTIME);
  return res;
}

