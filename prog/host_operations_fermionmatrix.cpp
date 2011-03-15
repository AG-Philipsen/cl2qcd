#include "host_operations_fermionmatrix.h"

hmc_error M(hmc_spinor_field* in, hmc_spinor_field* out, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im){
	M_diag(in, out, kappa, mu);    
	hmc_spinor_field tmp[SPINORFIELDSIZE];
	dslash(in,tmp,gaugefield, theta, chem_pot_re, chem_pot_im);

	hmc_complex kappa_cmplx = {kappa, 0.};
	saxpy(tmp, out, &kappa_cmplx, out);

	return HMC_SUCCESS;
}

hmc_error M_diag(hmc_spinor_field* in, hmc_spinor_field* out, hmc_float kappa, hmc_float mu){
	hmc_spinor spinout[SPINORSIZE];
	//iterate over all lattice sites
	for(int spacepos=0; spacepos<VOLSPACE; spacepos++) {
		for(int timepos=0; timepos<NTIME; timepos++) {
			get_spinor_from_field(in,spinout,spacepos,timepos);
			M_diag_local(spinout, kappa, mu);
			put_spinor_to_field(spinout,out,spacepos,timepos);
		}}
	return HMC_SUCCESS; 
}

hmc_error M_sitediagonal(hmc_eoprec_spinor_field* in, hmc_eoprec_spinor_field* out, hmc_float kappa, hmc_float mu){
	hmc_spinor spinout[SPINORSIZE];
	//iterate over half the lattice
  for(int n=0; n<VOL4D/2; n++) {
    get_spinor_from_eoprec_field(in,spinout,n);
    M_diag_local(spinout, kappa, mu);
    put_spinor_to_eoprec_field(spinout,out,n);
  }
  return HMC_SUCCESS;
}

hmc_error M_inverse_sitediagonal(hmc_eoprec_spinor_field* in, hmc_eoprec_spinor_field* out, hmc_float kappa, hmc_float mu){
	hmc_spinor spinout[SPINORSIZE];
	//iterate over half the lattice
	for(int n=0; n<VOL4D/2; n++) {
		hmc_float minuskappa = -kappa;
		get_spinor_from_eoprec_field(in,spinout,n);
		M_diag_local(spinout, minuskappa, mu);
		hmc_float denom = 1. + 4.*kappa*kappa*mu*mu;
		real_multiply_spinor(spinout,1./denom);
		put_spinor_to_eoprec_field(spinout,out,n);
	}
	return HMC_SUCCESS;
}

void dslash_spatial (hmc_spinor * spinout, int * coord, int dir, int pos, int t, hmc_spinor_field* in, hmc_gaugefield* gaugefield, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im){

	int next, prev;
	hmc_spinor spinnext[SPINORSIZE];
	hmc_spinor spinprev[SPINORSIZE];
	hmc_su3matrix u;
	hmc_su3matrix udagger; 

	next = get_neighbor(pos,dir);
	prev = get_lower_neighbor(pos,dir);
  
	get_spinor_from_field(in, spinnext, next, t);
	get_spinor_from_field(in, spinprev, prev, t);
      
	get_su3matrix(&u,gaugefield,pos,t,dir);
	get_su3matrix(&udagger,gaugefield,prev,t,dir);
	adjoin_su3matrix(&udagger);

	//update links with chemical potential, this shall be put into compiler option lateron
	gaugefield_apply_chem_pot(&u, &udagger, chem_pot_re, chem_pot_im);

	//!!CP: shouldnt this be only in time-direction??
	if(coord[dir] == NSPACE-1) spinor_apply_bc(spinnext, theta);
	else if(coord[dir] == 0) spinor_apply_bc(spinprev, theta);
      
	if(dir == 1) dslash_1(spinnext, spinprev, spinout, &u, &udagger);
	else if(dir == 2) dslash_2(spinnext, spinprev, spinout, &u, &udagger);
	else dslash_3(spinnext, spinprev, spinout, &u, &udagger);
	
	return;
}

void dslash_temporal (hmc_spinor * spinout, int pos, int t, hmc_spinor_field* in, hmc_gaugefield* gaugefield, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im){
	int next, prev;
	hmc_spinor spinnext[SPINORSIZE];
	hmc_spinor spinprev[SPINORSIZE];
	hmc_su3matrix u;
	hmc_su3matrix udagger; 

	next = (t+1)%NTIME; 
	prev = (t-1+NTIME)%NTIME;

	get_spinor_from_field(in, spinnext, pos, next);
	get_spinor_from_field(in, spinprev, pos, prev);

	if(next == 0) spinor_apply_bc(spinnext, theta);
	else if(prev == NTIME-1) spinor_apply_bc(spinprev, theta);
      
	get_su3matrix(&u,gaugefield,pos,t,0);
	get_su3matrix(&udagger,gaugefield,pos,prev,0);
	adjoin_su3matrix(&udagger);

	//update links with chemical potential, this shall be put into compiler option lateron
	gaugefield_apply_chem_pot(&u, &udagger, chem_pot_re, chem_pot_im);

	dslash_0(spinnext, spinprev, spinout, &u, &udagger);

	return;
}

void dslash_spatial_eoprec (hmc_spinor * spinout, int * coord, int dir, int pos, int t, hmc_eoprec_spinor_field* in, hmc_gaugefield* gaugefield, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im){

	int next, prev;
	int neo_next;
	int neo_prev;
	hmc_spinor spinnext[SPINORSIZE];
	hmc_spinor spinprev[SPINORSIZE];
	hmc_su3matrix u;
	hmc_su3matrix udagger; 

	next = get_neighbor(pos,dir);
	prev = get_lower_neighbor(pos,dir);

	neo_next = get_n_eoprec(t, next);
	neo_prev = get_n_eoprec(t, prev);

	get_spinor_from_eoprec_field(in,spinnext,neo_next);
	get_spinor_from_eoprec_field(in,spinprev,neo_prev);
      
	get_su3matrix(&u,gaugefield,pos,t,dir);
	get_su3matrix(&udagger,gaugefield,prev,t,dir);
	adjoin_su3matrix(&udagger);

	//update links with chemical potential, this shall be put into compiler option lateron
	gaugefield_apply_chem_pot(&u, &udagger, chem_pot_re, chem_pot_im);

	//!!CP: shouldnt this be only in time-direction??
	if(coord[dir] == NSPACE-1) spinor_apply_bc(spinnext, theta);
	else if(coord[dir] == 0) spinor_apply_bc(spinprev, theta);
      
	if(dir == 1) dslash_1(spinnext, spinprev, spinout, &u, &udagger);
	else if(dir == 2) dslash_2(spinnext, spinprev, spinout, &u, &udagger);
	else dslash_3(spinnext, spinprev, spinout, &u, &udagger);
	
	return;
}

void dslash_temporal_eoprec (hmc_spinor * spinout, int pos, int t, hmc_eoprec_spinor_field* in, hmc_gaugefield* gaugefield, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im){
	int next, prev;
	int neo_next;
	int neo_prev;
	hmc_spinor spinnext[SPINORSIZE];
	hmc_spinor spinprev[SPINORSIZE];
	hmc_su3matrix u;
	hmc_su3matrix udagger; 

	next = (t+1)%NTIME; 
	prev = (t-1+NTIME)%NTIME;

	neo_next = get_n_eoprec(pos, next);
	neo_prev = get_n_eoprec(pos, prev);
	
	get_spinor_from_eoprec_field(in,spinnext,neo_next);
	get_spinor_from_eoprec_field(in,spinprev,neo_prev);

	if(next == 0) spinor_apply_bc(spinnext, theta);
	else if(prev == NTIME-1) spinor_apply_bc(spinprev, theta);
      
	get_su3matrix(&u,gaugefield,pos,t,0);
	get_su3matrix(&udagger,gaugefield,pos,prev,0);
	adjoin_su3matrix(&udagger);

	//update links with chemical potential, this shall be put into compiler option lateron
	gaugefield_apply_chem_pot(&u, &udagger, chem_pot_re, chem_pot_im);

	dslash_0(spinnext, spinprev, spinout, &u, &udagger);

	return;
}

hmc_error dslash(hmc_spinor_field* in, hmc_spinor_field* out, hmc_gaugefield* gaugefield, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im){
	hmc_spinor spinout[SPINORSIZE];
	//iterate over the whole lattice
	for(int spacepos=0; spacepos<VOLSPACE; spacepos++) {
		for(int timepos=0; timepos<NTIME; timepos++) {
			set_local_zero_spinor(spinout);   
			//like in host_geometry
			int coord[NDIM];
			coord[0]=0;
			for(int j=1;j<NDIM;j++) coord[j] = get_spacecoord(spacepos,j);
			
			// spinout = U_0*(r-gamma_0)*spinnext + U^dagger_0(x-hat0) * (r+gamma_0)*spinprev
			dslash_temporal(spinout, spacepos, timepos, in, gaugefield, theta, chem_pot_re, chem_pot_im);
			// spinout += U_1*(r-gamma_1)*spinnext + U^dagger_1(x-hat1) * (r+gamma_1)*spinprev
			dslash_spatial (spinout, coord, 1, spacepos, timepos, in, gaugefield, theta, chem_pot_re, chem_pot_im);
			// spinout += U_2*(r-gamma_2)*spinnext + U^dagger_2(x-hat2) * (r+gamma_2)*spinprev
			dslash_spatial (spinout, coord, 2, spacepos, timepos, in, gaugefield, theta, chem_pot_re, chem_pot_im);
			// spinout += U_3*(r-gamma_3)*spinnext + U^dagger_3(x-hat3) * (r+gamma_3)*spinprev
			dslash_spatial (spinout, coord, 3, spacepos, timepos, in, gaugefield, theta, chem_pot_re, chem_pot_im);
			
			put_spinor_to_field(spinout,out,spacepos,timepos);
    }
  }

  return HMC_SUCCESS;
}

hmc_error dslash_eoprec(hmc_eoprec_spinor_field* in, hmc_eoprec_spinor_field* out, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im, int evenodd){
	hmc_spinor spinout[SPINORSIZE];
	int ns, nt;
	//iterate over half the lattice
  for(int n=0; n<VOL4D/2; n++) {
    set_local_zero_spinor(spinout);    

    if(evenodd == ODD) get_odd_site(n, &ns, &nt);
    else get_even_site(n, &ns, &nt);
       
    //like in host_geometry
    int coord[NDIM];
    coord[0]=0;
    for(int j=1;j<NDIM;j++) coord[j] = get_spacecoord(ns,j);    
    
		// spinout = U_0*(r-gamma_0)*spinnext + U^dagger_0(x-hat0) * (r+gamma_0)*spinprev
		dslash_temporal_eoprec (spinout, ns, nt, in, gaugefield, theta, chem_pot_re, chem_pot_im);
		// spinout += U_1*(r-gamma_1)*spinnext + U^dagger_1(x-hat1) * (r+gamma_1)*spinprev
		dslash_spatial_eoprec (spinout, coord, 1, ns, nt, in, gaugefield, theta, chem_pot_re, chem_pot_im);
		// spinout += U_2*(r-gamma_2)*spinnext + U^dagger_2(x-hat2) * (r+gamma_2)*spinprev
		dslash_spatial_eoprec (spinout, coord, 2, ns, nt, in, gaugefield, theta, chem_pot_re, chem_pot_im);   
		// spinout += U_3*(r-gamma_3)*spinnext + U^dagger_3(x-hat3) * (r+gamma_3)*spinprev
		dslash_spatial_eoprec (spinout, coord, 3, ns, nt, in, gaugefield, theta, chem_pot_re, chem_pot_im);
		
		//!!CP: why this kappa??
		real_multiply_spinor(spinout,-kappa);
		put_spinor_to_eoprec_field(spinout,out,n);
  }

  return HMC_SUCCESS;
}

//!!CP: the args have to be revisited because of change in the defs
hmc_error Aee(hmc_eoprec_spinor_field* in, hmc_eoprec_spinor_field* out, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float mu, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im){
  hmc_eoprec_spinor_field* spintmp1 = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field* spintmp2 = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field* spintmp3 = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field* spintmp4 = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_complex one = hmc_complex_one;
  set_zero_spinorfield_eoprec(spintmp3);
  
//CP: original 
  dslash_eoprec(in,spintmp1,gaugefield,kappa, theta, chem_pot_re, chem_pot_re,ODD); // D_oe
  M_inverse_sitediagonal(spintmp1,spintmp2,kappa,mu); // R_o^(-1)
  dslash_eoprec(spintmp2,out,gaugefield,kappa, theta, chem_pot_re, chem_pot_re,EVEN); // D_eo
  M_sitediagonal(in,spintmp1,kappa,mu); //R_e
  copy_spinor_eoprec(out, spintmp3);
  saxpy_eoprec(spintmp3, spintmp1, &one, out); 
  
  delete [] spintmp1;
  delete [] spintmp2;
  delete [] spintmp3;
  delete [] spintmp4;

  return HMC_SUCCESS;
}

/*
//CP: old M_diag
hmc_error M_diag(hmc_spinor_field* in, hmc_spinor_field* out, hmc_float kappa, hmc_float mu){
  for(int spacepos=0; spacepos<VOLSPACE; spacepos++) {
    for(int timepos=0; timepos<NTIME; timepos++) {
      hmc_spinor spinout[SPINORSIZE];
      hmc_spinor spintmp[SPINORSIZE];
      get_spinor_from_field(in,spinout,spacepos,timepos);
      hmc_float twistfactor = 2*kappa*mu;
      multiply_spinor_i_factor_gamma5(spinout,spintmp,twistfactor);
      spinors_accumulate(spinout,spintmp);
      put_spinor_to_field(spinout,out,spacepos,timepos);
    }
  }
  return HMC_SUCCESS; 
}

//CP: old M_sitediagonal
hmc_error M_sitediagonal(hmc_eoprec_spinor_field* out, hmc_eoprec_spinor_field* in, hmc_float kappa, hmc_float mu){
	hmc_spinor spinout[SPINORSIZE];
  for(int n=0; n<VOL4D/2; n++) {
    hmc_spinor tmp[SPINORSIZE];
    get_spinor_from_eoprec_field(in,spinout,n);
    hmc_float twistfactor = 2*kappa*mu;
    multiply_spinor_i_factor_gamma5(spinout,tmp,twistfactor);
    spinors_accumulate(spinout,tmp);
    put_spinor_to_eoprec_field(spinout,out,n);
  }
  return HMC_SUCCESS;
}

//CP: old M_inverse_sitediagonal
hmc_error M_inverse_sitediagonal(hmc_eoprec_spinor_field* out, hmc_eoprec_spinor_field* in, hmc_float kappa, hmc_float mu){
  for(int n=0; n<VOL4D/2; n++) {
    hmc_spinor tmp[SPINORSIZE];
    hmc_spinor spinout[SPINORSIZE];
    get_spinor_from_eoprec_field(in,spinout,n);
    hmc_float twistfactor = -2*kappa*mu;
    multiply_spinor_i_factor_gamma5(spinout,tmp,twistfactor);
    spinors_accumulate(spinout,tmp);
    hmc_float denom = 1 + twistfactor*twistfactor;
    real_multiply_spinor(spinout,1./denom);
    put_spinor_to_eoprec_field(spinout,out,n);
  }
  return HMC_SUCCESS;
}

//CP: old dslash
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

//CP: old dslash_eoprec
hmc_error dslash_eoprec(hmc_eoprec_spinor_field* out, hmc_eoprec_spinor_field* in, hmc_gaugefield* gaugefield, hmc_float kappa, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im, int evenodd){
	//iterate over half the lattice
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
     
    
//     spinprojectproduct_gamma0(&u,spinnext,-hmc_one_f);
//     spinors_accumulate(spinout,spinnext);
// 
//     spinprojectproduct_gamma0(&udagger,spinprev,hmc_one_f);
//     spinors_accumulate(spinout,spinprev);
    
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
      
    
//     spinprojectproduct_gamma1(&u,spinnext,-hmc_one_f);
//     spinors_accumulate(spinout,spinnext);
// 
//     spinprojectproduct_gamma1(&udagger,spinprev,hmc_one_f);
//     spinors_accumulate(spinout,spinprev);
    
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
    
    
//     spinprojectproduct_gamma2(&u,spinnext,-hmc_one_f);
//     spinors_accumulate(spinout,spinnext);
// 
//     spinprojectproduct_gamma2(&udagger,spinprev,hmc_one_f);
//     spinors_accumulate(spinout,spinprev);
    
    
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
    
    
//     spinprojectproduct_gamma3(&u,spinnext,-hmc_one_f);
//     spinors_accumulate(spinout,spinnext);
// 
//     spinprojectproduct_gamma3(&udagger,spinprev,hmc_one_f);
//     spinors_accumulate(spinout,spinprev);
     
    
    //CP: why this kappa??
    real_multiply_spinor(spinout,-kappa);
    put_spinor_to_eoprec_field(spinout,out,n);

  }

  return HMC_SUCCESS;
}

*/