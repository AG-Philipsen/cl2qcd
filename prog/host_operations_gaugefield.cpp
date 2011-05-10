#include "host_operations_gaugefield.h"

hmc_error copy_to_ocl_format(hmc_ocl_gaugefield* host_gaugefield,hmc_gaugefield* gaugefield){
  for(int spacepos=0; spacepos<NSPACE*NSPACE*NSPACE; spacepos++) {
    for(int t=0; t<NTIME; t++) {
      for(int mu=0; mu<NDIM; mu++) {
	hmc_su3matrix tmp;
	get_su3matrix(&tmp, gaugefield, spacepos, t, mu);
	for(int b=0; b<NC; b++) {
#ifdef _RECONSTRUCT_TWELVE_
	  for(int a=0; a<NC-1; a++) {
	    int n = a + (NC-1)*b;
	    host_gaugefield[ocl_gaugefield_element(0,a,b,mu,spacepos,t)] = tmp[n].re;
	    host_gaugefield[ocl_gaugefield_element(1,a,b,mu,spacepos,t)] = tmp[n].im;
	  }
#else
	  for(int a=0; a<NC; a++) {
	    host_gaugefield[ocl_gaugefield_element(0,a,b,mu,spacepos,t)] = tmp[a][b].re;
	    host_gaugefield[ocl_gaugefield_element(1,a,b,mu,spacepos,t)] = tmp[a][b].im;
	  }
#endif
	}
      }
    }
  }
  return HMC_SUCCESS;
}

hmc_error copy_from_ocl_format(hmc_gaugefield* gaugefield,hmc_ocl_gaugefield* host_gaugefield){
  for(int spacepos=0; spacepos<NSPACE*NSPACE*NSPACE; spacepos++) {
    for(int t=0; t<NTIME; t++) {
      for(int mu=0; mu<NDIM; mu++) {
				hmc_su3matrix tmp;
				for(int b=0; b<NC; b++) {
#ifdef _RECONSTRUCT_TWELVE_
	  			for(int a=0; a<NC-1; a++) {
	    			int n = a + (NC-1)*b;
	    			tmp[n].re = host_gaugefield[ocl_gaugefield_element(0,a,b,mu,spacepos,t)];
	    			tmp[n].im = host_gaugefield[ocl_gaugefield_element(1,a,b,mu,spacepos,t)];
	  			}
#else
	  			for(int a=0; a<NC; a++) {
						tmp[a][b].re = host_gaugefield[ocl_gaugefield_element(0,a,b,mu,spacepos,t)];
						tmp[a][b].im = host_gaugefield[ocl_gaugefield_element(1,a,b,mu,spacepos,t)];
	  			}
#endif
	  	put_su3matrix(gaugefield, &tmp, spacepos, t, mu);
				}}}}
  return HMC_SUCCESS;
}
 
hmc_error set_gaugefield_cold(hmc_gaugefield * field) {
  for(int t=0; t<NTIME; t++) {
    for(int n=0; n<VOLSPACE; n++) {
      for(int mu=0; mu<NDIM; mu++) {
	hmc_su3matrix tmp;
	unit_su3matrix(&tmp);
	put_su3matrix(field, &tmp, n, t, mu);
      }
    }
  }
  return HMC_SUCCESS;
}

hmc_error set_gaugefield_hot(hmc_gaugefield * field) {
  for(int t=0; t<NTIME; t++) {
    for(int n=0; n<VOLSPACE; n++) {
      for(int mu=0; mu<NDIM; mu++) {
	hmc_su3matrix tmp;
	random_su3matrix(&tmp);
	put_su3matrix(field, &tmp, n, t, mu);
      }
    }
  }
  return HMC_SUCCESS;
}

hmc_error copy_gaugefield_from_ildg_format(hmc_gaugefield * gaugefield, hmc_float * gaugefield_tmp, int check){
  //little check if arrays are big enough
  if (VOL4D*NDIM*NC*NC*2 != check){
    std::cout << "error in setting gaugefield to source values!! "<< std::endl << "Check global settings!!" << std::endl << std::endl;
    return HMC_STDERR;
  }

  int cter=0;
  //our def: hmc_gaugefield [NC][NC][NDIM][VOLSPACE][NTIME]([2]), last one implicit for complex
  for (int t = 0; t<NTIME; t++){
  //if the alg is known to be correct, the next three for-loops could also be changed to one
    for (int i = 0; i<NSPACE; i++){
      for (int j = 0; j<NSPACE; j++){
        for (int k = 0; k<NSPACE; k++){
          for (int l = 0; l<NDIM; l++){
            int spacepos = k + j*NSPACE + i*NSPACE*NSPACE;
            int globalpos = l + spacepos*NDIM + t*VOLSPACE*NDIM;
#ifdef _RECONSTRUCT_TWELVE_
            for (int m = 0; m<NC; m++){
              for (int n = 0; n<NC; n++){
		int ncindex = m + (NC-1)*n;
		//ildg-std: [NT][NZ][NY][NX][NDIMENSION][NCOLOR][NCOLOR][2]
		//which is stored in one single array here
		//skip NC*NC*2 cmplx numbers
		int pos = 2*n + 2*m*NC + globalpos*NC*NC*2;
		if(m<NC-1) {
		  (*gaugefield)[ncindex][(l+1)%NDIM][spacepos][t].re = gaugefield_tmp[pos];
		  (*gaugefield)[ncindex][(l+1)%NDIM][spacepos][t].im = gaugefield_tmp[pos + 1];
		}
		cter++;
	      }}
#else
            for (int m = 0; m<NC; m++){
              for (int n = 0; n<NC; n++){
                //ildg-std: [NT][NZ][NY][NX][NDIMENSION][NCOLOR][NCOLOR][2]
                //which is stored in one single array here
                //skip NC*NC*2 cmplx numbers
                int pos = 2*n + 2*m*NC + globalpos*NC*NC*2;
                (*gaugefield)[m][n][(l+1)%NDIM][spacepos][t].re = gaugefield_tmp[pos];
                (*gaugefield)[m][n][(l+1)%NDIM][spacepos][t].im = gaugefield_tmp[pos + 1];
                cter++;
	      }}
#endif
  }}}}}

  if(cter*2 != check) {
    std::cout << "error in setting gaugefield to source values! there were " << cter*2 << " vals set and not " << check << std::endl;
    return HMC_STDERR;
  }

  return HMC_SUCCESS;
}

hmc_error copy_gaugefield_to_ildg_format(ildg_gaugefield * dest, hmc_gaugefield * source){
  
  int cter=0;
  //our def: hmc_gaugefield [NC][NC][NDIM][VOLSPACE][NTIME]([2]), last one implicit for complex
  for (int t = 0; t<NTIME; t++){
  //if the alg is known to be correct, the next three for-loops could also be changed to one
    for (int i = 0; i<NSPACE; i++){
      for (int j = 0; j<NSPACE; j++){
        for (int k = 0; k<NSPACE; k++){
          for (int l = 0; l<NDIM; l++){
            int spacepos = k + j*NSPACE + i*NSPACE*NSPACE;
            int globalpos = l + spacepos*NDIM + t*VOLSPACE*NDIM;
#ifdef _RECONSTRUCT_TWELVE_
            for (int m = 0; m<NC; m++){
              for (int n = 0; n<NC; n++){
		int ncindex = m + (NC-1)*n;
		//ildg-std: [NT][NZ][NY][NX][NDIMENSION][NCOLOR][NCOLOR][2]
		//which is stored in one single array here
		//skip NC*NC*2 cmplx numbers
		int pos = 2*n + 2*m*NC + globalpos*NC*NC*2;
		if(m<NC-1) {
		  (dest[0])[pos]     = (*source)[ncindex][(l+1)%NDIM][spacepos][t].re;
		  (dest[0])[pos + 1] = (*source)[ncindex][(l+1)%NDIM][spacepos][t].im;
		}
		cter++;
	      }}
#else
            for (int m = 0; m<NC; m++){
              for (int n = 0; n<NC; n++){
                //ildg-std: [NT][NZ][NY][NX][NDIMENSION][NCOLOR][NCOLOR][2]
                //which is stored in one single array here
                //skip NC*NC*2 cmplx numbers
                int pos = 2*n + 2*m*NC + globalpos*NC*NC*2;
                (dest[0])[pos]     = ((*source)[m][n][(l+1)%NDIM][spacepos][t]).re;
                (dest[0])[pos + 1] = ((*source)[m][n][(l+1)%NDIM][spacepos][t]).im;
                cter++;
	      }}
#endif
  }}}}}
  
  return HMC_SUCCESS;
}

hmc_complex global_trace_su3(hmc_gaugefield * field, int mu) {
  hmc_complex sum;
  sum.re = 0;
  sum.im = 0;
  for(int t=0; t<NTIME; t++) {
    for(int n=0; n<VOLSPACE; n++) {
      hmc_su3matrix tmp;
      get_su3matrix(&tmp, field, n, t, mu);
      sum.re += trace_su3matrix(&tmp).re;
      sum.im += trace_su3matrix(&tmp).im;
    }
  }
  return sum;
}




hmc_error get_su3matrix(hmc_su3matrix * out, hmc_gaugefield * in, int spacepos, int timepos, int mu) {
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*(NC-1); n++) (*out)[n] = (*in)[n][mu][spacepos][timepos];
#else
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      (*out)[a][b] = (*in)[a][b][mu][spacepos][timepos];
    }
  }
#endif
  return HMC_SUCCESS;
}

hmc_error put_su3matrix(hmc_gaugefield * field, hmc_su3matrix * in, int spacepos, int timepos, int mu) {
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*(NC-1); n++) (*field)[n][mu][spacepos][timepos] = (*in)[n];
#else
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      (*field)[a][b][mu][spacepos][timepos] = (*in)[a][b];
    }
  }
#endif
  return HMC_SUCCESS;
}

void local_polyakov(hmc_gaugefield * field, hmc_su3matrix * prod, int n){
	unit_su3matrix(prod);
	for(int t=0; t<NTIME; t++) {
		hmc_su3matrix tmp;
		get_su3matrix(&tmp,field,n,t,0);
		accumulate_su3matrix_prod(prod,&tmp);
	}
	return;
}

void local_plaquette(hmc_gaugefield * field, hmc_su3matrix * prod, int n, int t, int mu, int nu ){
	hmc_su3matrix tmp;
	//u_mu(x)
	get_su3matrix(prod,field,n,t,mu);
	//u_nu(x+mu)
	if(mu==0) {
	  int newt = (t+1)%NTIME; //(haha)
	  get_su3matrix(&tmp,field,n,newt,nu);
	} else {
	  get_su3matrix(&tmp,field,get_neighbor(n,mu),t,nu);
	}
	accumulate_su3matrix_prod(prod,&tmp);
	//adjoint(u_mu(x+nu))
	if(nu==0) {
	  int newt = (t+1)%NTIME;
	  get_su3matrix(&tmp,field,n,newt,mu);
	} else {
	  get_su3matrix(&tmp,field,get_neighbor(n,nu),t,mu);
	}
	adjoin_su3matrix(&tmp);
	accumulate_su3matrix_prod(prod,&tmp);
	//adjoint(u_nu(x))
	get_su3matrix(&tmp,field,n,t,nu);
	adjoin_su3matrix(&tmp);
	accumulate_su3matrix_prod(prod,&tmp);
		
	return;
}

/** @todo memcpy ... */
hmc_error copy_gaugefield(hmc_gaugefield * source, hmc_gaugefield * dest){
	// copies source to destination within cpu memory, layer for gaugefield array
	return complexcopy((hmc_complex *)source, (hmc_complex *)dest, GAUGEFIELDSIZE); // SL: not tested
}

