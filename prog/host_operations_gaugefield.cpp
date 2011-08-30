#include "host_operations_gaugefield.h"

void copy_to_ocl_format(ocl_s_gaugefield* host_gaugefield, s_gaugefield* gaugefield){
  for(int spacepos=0; spacepos<NSPACE*NSPACE*NSPACE; spacepos++) {
    for(int t=0; t<NTIME; t++) {
      for(int mu=0; mu<NDIM; mu++) {
	host_gaugefield[get_global_link_pos(mu, spacepos, t)] = (*gaugefield) [mu][spacepos][t];
      }
    }
  }
  return;
}

void copy_from_ocl_format(s_gaugefield* gaugefield, ocl_s_gaugefield* host_gaugefield){
  for(int spacepos=0; spacepos<NSPACE*NSPACE*NSPACE; spacepos++) {
    for(int t=0; t<NTIME; t++) {
      for(int mu=0; mu<NDIM; mu++) {
 	(*gaugefield) [mu][spacepos][t] = host_gaugefield [get_global_link_pos(mu, spacepos, t)];
      }
    }
  }
  return;
}

void set_gaugefield_cold(hmc_gaugefield * field) {
  for(int t=0; t<NTIME; t++) {
    for(int n=0; n<VOLSPACE; n++) {
      for(int mu=0; mu<NDIM; mu++) {
	hmc_su3matrix tmp;
	unit_su3matrix(&tmp);
	put_su3matrix(field, &tmp, n, t, mu);
      }
    }
  }
  return;
}

void set_gaugefield_hot(hmc_gaugefield * field) {
  for(int t=0; t<NTIME; t++) {
    for(int n=0; n<VOLSPACE; n++) {
      for(int mu=0; mu<NDIM; mu++) {
	hmc_su3matrix tmp;
	random_su3matrix(&tmp);
	put_su3matrix(field, &tmp, n, t, mu);
      }
    }
  }
  return;
}

void copy_gaugefield_from_ildg_format(hmc_gaugefield * gaugefield, hmc_float * gaugefield_tmp, int check){
  //little check if arrays are big enough
  if (VOL4D*NDIM*NC*NC*2 != check){
    std::stringstream errstr;
    errstr<<"Error in setting gaugefield to source values!!\nCheck global settings!!";
    throw Print_Error_Message(errstr.str(),__FILE__,__LINE__);
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
    std::stringstream errstr;
    errstr << "Error in setting gaugefield to source values! there were " << cter*2 << " vals set and not " << check <<".";
    throw Print_Error_Message(errstr.str(),__FILE__,__LINE__);
  }

  return;
}

void copy_gaugefield_to_ildg_format(ildg_gaugefield * dest, hmc_gaugefield * source){
  
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
		if (m==NC-1){
		  hmc_su3matrix tmpsu3;
		  get_su3matrix(&tmpsu3, source, spacepos, t, (l+1)%NDIM);
		  hmc_complex tmp = reconstruct_su3 (&tmpsu3, n);
		  (dest[0])[pos]     = tmp.re;
		  (dest[0])[pos + 1] = tmp.im;
		}
		
// 		if ((dest[0])[pos] != (dest[0])[pos])
// 		  cout << (dest[0])[pos] <<endl;
// 		if ((dest[0])[pos+1] != (dest[0])[pos+1])
// 		  cout << (dest[0])[pos+1]   <<endl;
		
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
  
  return;
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




void get_su3matrix(hmc_su3matrix * out, hmc_gaugefield * in, int spacepos, int timepos, int mu) {
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*(NC-1); n++) (*out)[n] = (*in)[n][mu][spacepos][timepos];
#else
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      (*out)[a][b] = (*in)[a][b][mu][spacepos][timepos];
    }
  }
#endif
  return;
}

void put_su3matrix(hmc_gaugefield * field, hmc_su3matrix * in, int spacepos, int timepos, int mu) {
#ifdef _RECONSTRUCT_TWELVE_
  for(int n=0; n<NC*(NC-1); n++) (*field)[n][mu][spacepos][timepos] = (*in)[n];
#else
  for(int a=0; a<NC; a++) {
    for(int b=0; b<NC; b++) {
      (*field)[a][b][mu][spacepos][timepos] = (*in)[a][b];
    }
  }
#endif
  return;
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
void copy_gaugefield(hmc_gaugefield * source, hmc_gaugefield * dest){
	// copies source to destination within cpu memory, layer for gaugefield array
	return complexcopy((hmc_complex *)source, (hmc_complex *)dest, GAUGEFIELDSIZE); // SL: not tested
}

void local_Q_plaquette(hmc_3x3matrix * out, hmc_gaugefield * field, int n, int t, int mu, int nu ){
  hmc_su3matrix tmp;
  int newpos;
  
  //1st plaq
  hmc_su3matrix plaq1;
  //u_mu(x)
  get_su3matrix(&plaq1,field,n,t,mu);
  //u_nu(x+mu)
  if(mu==0) {
    int newt = (t+1)%NTIME;
    get_su3matrix(&tmp,field,n,newt,nu);
  }
  else
    get_su3matrix(&tmp,field,get_neighbor(n,mu),t,nu);
  accumulate_su3matrix_prod(&plaq1,&tmp);
  //adjoint(u_mu(x+nu))
  if(nu==0) {
    int newt = (t+1)%NTIME;
    get_su3matrix(&tmp,field,n,newt,mu);
  }
  else
    get_su3matrix(&tmp,field,get_neighbor(n,nu),t,mu);
  adjoin_su3matrix(&tmp);
  accumulate_su3matrix_prod(&plaq1,&tmp);
  //adjoint(u_nu(x))
  get_su3matrix(&tmp,field,n,t,nu);
  adjoin_su3matrix(&tmp);
  accumulate_su3matrix_prod(&plaq1,&tmp);

  //2nd plaq
  hmc_su3matrix plaq2;
  //U_nu(x)
  get_su3matrix(&plaq2,field,n,t,nu);
  //adj (u_mu(x-mu+nu))
  newpos = get_lower_neighbor(n, mu);
  if (nu==0){
    int newt =  (t+1)%NTIME;
    get_su3matrix(&tmp,field,newpos,newt,mu);
  }
  else if (mu==0){
    int newt = (t-1+NTIME)%NTIME;
    get_su3matrix(&tmp,field,get_neighbor(n,nu),newt,mu);
  }
  else
    get_su3matrix(&tmp,field,get_neighbor(newpos, nu),t,mu);
  adjoin_su3matrix(&tmp);
  accumulate_su3matrix_prod(&plaq2,&tmp);
  //adj (u_nu(x-mu))
  if (mu==0){
    int newt = (t-1+NTIME)%NTIME;
    get_su3matrix(&tmp,field,n,newt, nu);
  }
  else
    get_su3matrix(&tmp,field, get_lower_neighbor(n, mu),t, nu);
  adjoin_su3matrix(&tmp);
  accumulate_su3matrix_prod(&plaq2,&tmp);
  //u_mu(x-mu)
  if (mu==0){
    int newt = (t-1+NTIME)%NTIME;
    get_su3matrix(&tmp,field,n,newt, mu);
  }
  else
    get_su3matrix(&tmp,field, get_lower_neighbor(n, mu),t, mu);
  accumulate_su3matrix_prod(&plaq2,&tmp);

  //3rd plaq
  hmc_su3matrix plaq3;
  //adj (u_mu(x-mu))
  if (mu==0){
    int newt = (t-1+NTIME)%NTIME;
    get_su3matrix(&tmp,field,n,newt, mu);
  }
  else
    get_su3matrix(&tmp,field, get_lower_neighbor(n, mu),t, mu);
  adjoin_su3matrix(&tmp);
  copy_su3matrix(&plaq3, &tmp); 
  //adj (u_nu(x-mu-nu))
  newpos = get_lower_neighbor(n, mu);
  if (nu==0){
    int newt = (t-1+NTIME)%NTIME;
    get_su3matrix(&tmp,field, newpos,newt,nu);
  }
  else if (mu==0){
    int newt = (t-1+NTIME)%NTIME;
  get_su3matrix(&tmp,field,get_lower_neighbor(n,nu),newt,nu);
  }
  else
    get_su3matrix(&tmp,field,get_lower_neighbor(newpos, nu),t,nu);
  adjoin_su3matrix(&tmp);
  accumulate_su3matrix_prod(&plaq3,&tmp);
  //u_mu(x-mu-nu)
  newpos = get_lower_neighbor(n, mu);
  if (nu==0){
    int newt = (t-1+NTIME)%NTIME;
    get_su3matrix(&tmp,field,newpos,newt,mu);
  }
  else if (mu==0){
    int newt = (t-1+NTIME)%NTIME;
    get_su3matrix(&tmp,field,get_lower_neighbor(n,nu),newt,mu);
  }
  else
    get_su3matrix(&tmp,field,get_lower_neighbor(newpos, nu),t,mu);
  accumulate_su3matrix_prod(&plaq3,&tmp);
  //u_nu(x-nu)
  if (nu==0){
    int newt = (t-1+NTIME)%NTIME;
    get_su3matrix(&tmp,field, n,newt, nu);
  }
  else
    get_su3matrix(&tmp,field,get_lower_neighbor(n, nu),t,nu);
  accumulate_su3matrix_prod(&plaq3,&tmp);

  //4th plaq
  hmc_su3matrix plaq4;
  //adj(u_nu(x-nu))
  adjoin_su3matrix(&tmp);
  copy_su3matrix(&plaq4, &tmp); 
  //u_mu(x-nu)
   if (nu==0){
    int newt = (t-1+NTIME)%NTIME;
    get_su3matrix(&tmp,field, n,newt, mu);
  }
  else
    get_su3matrix(&tmp,field, get_lower_neighbor(n, nu),t,mu);
  accumulate_su3matrix_prod(&plaq4,&tmp);
  //u_nu(x+mu-nu)
  newpos = get_lower_neighbor(n, nu);
  if (mu==0){
    int newt =  (t+1)%NTIME;
    get_su3matrix(&tmp,field,newpos,newt,nu);
  }
  else if (nu==0){
    int newt = (t-1+NTIME)%NTIME;
    get_su3matrix(&tmp,field,get_neighbor(n,mu),newt,nu);
  }
  else
    get_su3matrix(&tmp,field,get_neighbor(newpos, mu),t,nu);
  accumulate_su3matrix_prod(&plaq4,&tmp);
  //adj (u_mu(x))
  get_su3matrix(&tmp,field,n,t,mu);
  adjoin_su3matrix(&tmp);
  accumulate_su3matrix_prod(&plaq4,&tmp);

  //Sum up
  su3matrix_to_3x3matrix (out, &plaq1);
  accumulate_su3matrix_3x3_add(out, &plaq2);
  accumulate_su3matrix_3x3_add(out, &plaq3);
  accumulate_su3matrix_3x3_add(out, &plaq4);
}
