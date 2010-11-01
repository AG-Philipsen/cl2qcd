#include "gaugeobservables.h"
#include <cstdio>

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
	    int newt = (t+mu)%NTIME; //(haha)
	    get_su3matrix(&tmp,field,n,newt,mu);
	  } else {
	    get_su3matrix(&tmp,field,get_neighbor(n,mu),t,mu);
	  }
	  accumulate_su3matrix_prod(&prod,&tmp);
	  //adjoint(u_mu(x+nu))
	  if(nu==0) {
	    int newt = (t+nu)%NTIME;
	    get_su3matrix(&tmp,field,n,newt,nu);
	  } else {
	    get_su3matrix(&tmp,field,get_neighbor(n,nu),t,nu);
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
  *splaq /= 2.0 / static_cast<hmc_float>(VOL4D*NC*(NDIM-1)*(NDIM-2));
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
