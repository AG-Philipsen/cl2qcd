#include "host_operations_fermionmatrix.h"

hmc_error M(inputparameters * parameters, hmc_spinor_field* in, hmc_gaugefield* gaugefield, hmc_spinor_field* out){
	M_diag(parameters, in, out);    
	hmc_spinor_field tmp[SPINORFIELDSIZE];
	dslash(parameters, in,gaugefield,tmp);

	hmc_complex kappa_cmplx = {(*parameters).get_kappa(), 0.};
	saxpy(tmp, out, &kappa_cmplx, out);

	return HMC_SUCCESS;
}

hmc_error Mdagger(inputparameters * parameters, hmc_spinor_field* in, hmc_gaugefield* gaugefield, hmc_spinor_field* out){
	Mdagger_diag(parameters, in, out);    
	hmc_spinor_field tmp[SPINORFIELDSIZE];
	ddaggerslash(parameters, in,gaugefield,tmp);

	hmc_complex kappa_cmplx = {(*parameters).get_kappa(), 0.};
	saxpy(tmp, out, &kappa_cmplx, out);

	return HMC_SUCCESS;
}

hmc_error M_diag(inputparameters * parameters, hmc_spinor_field* in, hmc_spinor_field* out){
	hmc_spinor spinout[SPINORSIZE];
	//iterate over all lattice sites
	for(int spacepos=0; spacepos<VOLSPACE; spacepos++) {
		for(int timepos=0; timepos<NTIME; timepos++) {
			get_spinor_from_field(in,spinout,spacepos,timepos);
			M_diag_local(spinout, (*parameters).get_kappa(), (*parameters).get_mu());
			put_spinor_to_field(spinout,out,spacepos,timepos);
		}}
	return HMC_SUCCESS; 
}

//if one would have explicit flavour structure, this would have to be revisited because (pauli)dagger is non-trivial
hmc_error Mdagger_diag(inputparameters * parameters, hmc_spinor_field* in, hmc_spinor_field* out){
	hmc_spinor spinout[SPINORSIZE];
	//iterate over all lattice sites
	for(int spacepos=0; spacepos<VOLSPACE; spacepos++) {
		for(int timepos=0; timepos<NTIME; timepos++) {
			get_spinor_from_field(in,spinout,spacepos,timepos);
			hmc_float tmp = -(*parameters).get_kappa();
			M_diag_local(spinout, tmp, (*parameters).get_mu());
			put_spinor_to_field(spinout,out,spacepos,timepos);
		}}
	return HMC_SUCCESS; 
}

hmc_error MdaggerM_diag(inputparameters * parameters, hmc_spinor_field* in, hmc_spinor_field* out){
	hmc_spinor spinout[SPINORSIZE];
	//iterate over all lattice sites
	for(int spacepos=0; spacepos<VOLSPACE; spacepos++) {
		for(int timepos=0; timepos<NTIME; timepos++) {
			get_spinor_from_field(in,spinout,spacepos,timepos);
			MdaggerM_diag_local(spinout, (*parameters).get_kappa(), (*parameters).get_mu());
			put_spinor_to_field(spinout,out,spacepos,timepos);
		}}
	return HMC_SUCCESS; 
}

hmc_error MdaggerM(inputparameters * parameters, hmc_spinor_field* in, hmc_gaugefield* gaugefield, hmc_spinor_field* out){
	hmc_spinor_field tmp[SPINORFIELDSIZE];
	hmc_spinor_field tmp2[SPINORFIELDSIZE];
	hmc_complex kappa_cmplx = {(*parameters).get_kappa(), 0.};
	
	//diagonal part: out = MdaggerM_diag
	MdaggerM_diag(parameters, in, out);    
	//out += - kappa*Ddagger( Mdiag in )
	M_diag(parameters, in, tmp);
	ddaggerslash(parameters, tmp,gaugefield,tmp2);
	saxpy(tmp2, out, &kappa_cmplx, out);
	//out += - kappa*D( Mdaggerdiag in )
	Mdagger_diag(parameters, in, tmp);
	dslash(parameters, tmp,gaugefield,tmp2);
	saxpy(tmp2, out, &kappa_cmplx, out);
	//out += kappa*kappa* DdaggerD
	kappa_cmplx = {-(*parameters).get_kappa()*(*parameters).get_kappa(),0.};
	ddaggerd(parameters, in,gaugefield,tmp2);
	saxpy(tmp2, out, &kappa_cmplx, out);
	return HMC_SUCCESS;
}

hmc_error M_sitediagonal(inputparameters * parameters, hmc_eoprec_spinor_field* in, hmc_eoprec_spinor_field* out){
	hmc_spinor spinout[SPINORSIZE];
	//iterate over half the lattice
  for(int n=0; n<VOL4D/2; n++) {
    get_spinor_from_eoprec_field(in,spinout,n);
    M_diag_local(spinout, (*parameters).get_kappa(), (*parameters).get_mu());
    put_spinor_to_eoprec_field(spinout,out,n);
  }
  return HMC_SUCCESS;
}

hmc_error M_inverse_sitediagonal(inputparameters * parameters, hmc_eoprec_spinor_field* in, hmc_eoprec_spinor_field* out){
	
	hmc_spinor spinout[SPINORSIZE];
	//iterate over half the lattice
	for(int n=0; n<VOL4D/2; n++) {
		hmc_float minuskappa = -(*parameters).get_kappa();
		get_spinor_from_eoprec_field(in,spinout,n);
		M_diag_local(spinout, minuskappa, (*parameters).get_mu());
		hmc_float denom = 1. + 4.*(*parameters).get_kappa()*(*parameters).get_kappa()*(*parameters).get_mu()*(*parameters).get_mu();
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

	// /todo CP: shouldnt this be only in time-direction??
	if(coord[dir] == NSPACE-1) spinor_apply_bc(spinnext, theta);
	else if(coord[dir] == 0) spinor_apply_bc(spinprev, theta);
      
	if(dir == 1) dslash_1(spinnext, spinprev, spinout, &u, &udagger);
	else if(dir == 2) dslash_2(spinnext, spinprev, spinout, &u, &udagger);
	else dslash_3(spinnext, spinprev, spinout, &u, &udagger);
	
	return;
}

void ddaggerslash_spatial (hmc_spinor * spinout, int * coord, int dir, int pos, int t, hmc_spinor_field* in, hmc_gaugefield* gaugefield, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im){

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
      
	if(dir == 1) ddaggerslash_1(spinnext, spinprev, spinout, &u, &udagger);
	else if(dir == 2) ddaggerslash_2(spinnext, spinprev, spinout, &u, &udagger);
	else ddaggerslash_3(spinnext, spinprev, spinout, &u, &udagger);
	
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

void ddaggerslash_temporal (hmc_spinor * spinout, int pos, int t, hmc_spinor_field* in, hmc_gaugefield* gaugefield, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im){
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

	ddaggerslash_0(spinnext, spinprev, spinout, &u, &udagger);

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

	neo_next = get_n_eoprec(next, t);
	neo_prev = get_n_eoprec(prev, t);

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

hmc_error dslash(inputparameters * parameters, hmc_spinor_field* in, hmc_gaugefield* gaugefield, hmc_spinor_field* out){
	hmc_float kappa; hmc_float mu; hmc_float theta; hmc_float chem_pot_re; hmc_float chem_pot_im; int cgmax;
  kappa = (*parameters).get_kappa(); mu = (*parameters).get_mu(); theta = (*parameters).get_theta_fermion(); chem_pot_re = (*parameters).get_chem_pot_re(); chem_pot_im = (*parameters).get_chem_pot_im(); cgmax = (*parameters).get_cgmax();
		
	
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

hmc_error ddaggerslash(inputparameters * parameters, hmc_spinor_field* in, hmc_gaugefield* gaugefield, hmc_spinor_field* out){
	hmc_float kappa; hmc_float mu; hmc_float theta; hmc_float chem_pot_re; hmc_float chem_pot_im; int cgmax;
  kappa = (*parameters).get_kappa(); mu = (*parameters).get_mu(); theta = (*parameters).get_theta_fermion(); chem_pot_re = (*parameters).get_chem_pot_re(); chem_pot_im = (*parameters).get_chem_pot_im(); cgmax = (*parameters).get_cgmax();
		
	hmc_spinor spinout[SPINORSIZE];
	//iterate over the whole lattice
	for(int spacepos=0; spacepos<VOLSPACE; spacepos++) {
		for(int timepos=0; timepos<NTIME; timepos++) {
			set_local_zero_spinor(spinout);   
			//like in host_geometry
			int coord[NDIM];
			coord[0]=0;
			for(int j=1;j<NDIM;j++) coord[j] = get_spacecoord(spacepos,j);
			
			// spinout = U_0*(r+gamma_0)*spinnext + U^dagger_0(x-hat0) * (r-gamma_0)*spinprev
			ddaggerslash_temporal(spinout, spacepos, timepos, in, gaugefield, theta, chem_pot_re, chem_pot_im);
			// spinout += U_1*(r+gamma_1)*spinnext + U^dagger_1(x-hat1) * (r-gamma_1)*spinprev
			ddaggerslash_spatial (spinout, coord, 1, spacepos, timepos, in, gaugefield, theta, chem_pot_re, chem_pot_im);
			// spinout += U_2*(r+gamma_2)*spinnext + U^dagger_2(x-hat2) * (r-gamma_2)*spinprev
			ddaggerslash_spatial (spinout, coord, 2, spacepos, timepos, in, gaugefield, theta, chem_pot_re, chem_pot_im);
			// spinout += U_3*(r+gamma_3)*spinnext + U^dagger_3(x-hat3) * (r-gamma_3)*spinprev
			ddaggerslash_spatial (spinout, coord, 3, spacepos, timepos, in, gaugefield, theta, chem_pot_re, chem_pot_im);
			
			put_spinor_to_field(spinout,out,spacepos,timepos);
    }
  }

  return HMC_SUCCESS;
}

hmc_error ddaggerd(inputparameters * parameters, hmc_spinor_field* in, hmc_gaugefield* gaugefield, hmc_spinor_field* out){
	hmc_float kappa; hmc_float mu; hmc_float theta; hmc_float chem_pot_re; hmc_float chem_pot_im; int cgmax;
  kappa = (*parameters).get_kappa(); mu = (*parameters).get_mu(); theta = (*parameters).get_theta_fermion(); chem_pot_re = (*parameters).get_chem_pot_re(); chem_pot_im = (*parameters).get_chem_pot_im(); cgmax = (*parameters).get_cgmax();
		
	hmc_spinor spinout[SPINORSIZE];
	//iterate over the whole lattice
	for(int spacepos=0; spacepos<VOLSPACE; spacepos++) {
		for(int timepos=0; timepos<NTIME; timepos++) {
			set_local_zero_spinor(spinout);   
			
			ddaggerd_calc(spinout, spacepos, timepos, in, gaugefield, theta, chem_pot_re, chem_pot_im);
		
			put_spinor_to_field(spinout,out,spacepos,timepos);
    }
  }

  return HMC_SUCCESS;
}

hmc_error dslash_eoprec(inputparameters * parameters, hmc_eoprec_spinor_field* in, hmc_gaugefield* gaugefield, int evenodd, hmc_eoprec_spinor_field* out){
	hmc_float kappa; hmc_float mu; hmc_float theta; hmc_float chem_pot_re; hmc_float chem_pot_im; int cgmax;
  kappa = (*parameters).get_kappa(); mu = (*parameters).get_mu(); theta = (*parameters).get_theta_fermion(); chem_pot_re = (*parameters).get_chem_pot_re(); chem_pot_im = (*parameters).get_chem_pot_im(); cgmax = (*parameters).get_cgmax();	
	
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

hmc_error Aee(inputparameters * parameters, hmc_eoprec_spinor_field* in, hmc_gaugefield* gaugefield, hmc_eoprec_spinor_field* out){

	hmc_eoprec_spinor_field* spintmp1 = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field* spintmp2 = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_eoprec_spinor_field* spintmp3 = new hmc_eoprec_spinor_field[EOPREC_SPINORFIELDSIZE];
  hmc_complex one = hmc_complex_one;
  set_zero_spinorfield_eoprec(spintmp3);
  
  dslash_eoprec(parameters, in,gaugefield,ODD, spintmp1); // D_oe
  M_inverse_sitediagonal(parameters, spintmp1,spintmp2); // R_o^(-1)
  dslash_eoprec(parameters, spintmp2,gaugefield,EVEN, out); // D_eo
  M_sitediagonal(parameters, in,spintmp1); //R_e
  copy_spinor_eoprec(out, spintmp3);
  saxpy_eoprec(spintmp3, spintmp1, &one, out); 
  
  delete [] spintmp1;
  delete [] spintmp2;
  delete [] spintmp3;

  return HMC_SUCCESS;
}


void ddaggerd_calc (hmc_spinor * spinout, int pos, int t, hmc_spinor_field* in, hmc_gaugefield* gaugefield, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im){
	int next_time, next_spat, prev_time, prev_spat;
	int nextplusnu_time, nextminusnu_time;
	int prevplusnu_time, prevminusnu_time;
	int nextplusnu_spat, nextminusnu_spat;
	int prevplusnu_spat, prevminusnu_spat;
	
	hmc_spinor spinnext[SPINORSIZE];
	hmc_spinor spinprev[SPINORSIZE];
	hmc_spinor tmp1[SPINORSIZE];
	hmc_spinor tmp2[SPINORSIZE];
	hmc_su3matrix u;
	hmc_su3matrix udagger; 
	
	int nu[3], nu_spat[3], nu_time[3];

	int mu;
	//outer loop
	for (mu = 0; mu<NDIM; mu++)
	{

#define SET3( dest, val1, val2, val3 ) \
{ \
	dest[0] = val1; \
	dest[1] = val2; \
	dest[2] = val3; \
}

		//directions where calculations have to be carried out
		if(mu == 0){
			SET3(nu, 1,2,3);
			SET3(nu_spat, 1,2,3);
			SET3(nu_time, 0,0,0);
		}
		else if (mu == 1){
			SET3(nu, 0,2,3);
			SET3(nu_spat, 0,2,3);
			SET3(nu_time, 1,0,0);
		}
		else if (mu == 2){
			SET3(nu, 0,1,3);
			SET3(nu_spat, 0,1,3);
			SET3(nu_time, 1,0,0);
		}
		else{
			SET3(nu, 0,1,2);
			SET3(nu_spat, 0,1,2);
			SET3(nu_time, 1,0,0);
		}

#undef SET3

		//get new positions
		if(mu == 0){
			next_time = (t+1)%NTIME; 
			prev_time = (t-1+NTIME)%NTIME;
			next_spat = pos;
			prev_spat = pos;
		}
		else{
			next_time = t; 
			prev_time = t;
			next_spat = get_neighbor(pos, mu);
			prev_spat = get_lower_neighbor(pos, mu);
		}
		
		set_local_zero_spinor(tmp1);
		set_local_zero_spinor(tmp2);
		//iterate through nu neq mu
		for(int i = 0; i<NDIM-1; i++){
			//calcl coordinates needed
			//n + mu + nu
			nextplusnu_spat = get_neighbor(next_spat,nu_spat[i]);
			nextplusnu_time = (next_time + nu_time[i] )%NTIME;
			//n + mu - nu
			nextminusnu_spat = get_lower_neighbor(next_spat,nu_spat[i]);
			nextplusnu_time = (next_time - nu_time[i] )%NTIME;
			//n - mu + nu
			prevplusnu_spat = get_neighbor(prev_spat,nu_spat[i]);
			prevplusnu_time = (prev_time + nu_time[i] )%NTIME;
			//n - mu - nu
			prevminusnu_spat = get_lower_neighbor(prev_spat,nu_spat[i]);
			prevminusnu_time = (prev_time - nu_time[i] )%NTIME;
		
			//first two
			get_spinor_from_field(in, spinnext, prevplusnu_spat, prevplusnu_time);
			get_spinor_from_field(in, spinprev, prevminusnu_spat, prevminusnu_time);

			//TODO
// 		if(prev_time == 0) spinor_apply_bc(spinnext, theta);
// 		else if(prevminusnu_time == NTIME-1) spinor_apply_bc(spinprev, theta);
      
			get_su3matrix(&u,gaugefield,prev_spat, prev_time,nu[i]);
			get_su3matrix(&udagger,gaugefield,prevminusnu_spat,prevminusnu_time,mu);
			adjoin_su3matrix(&udagger);

			//update links with chemical potential, this shall be put into compiler option lateron
			gaugefield_apply_chem_pot(&u, &udagger, chem_pot_re, chem_pot_im);

			//store in tmp1
			if(mu == 0){
				dslash_0(spinnext, spinprev, tmp1, &u, &udagger);
			}
			else if (mu == 1){
				dslash_1(spinnext, spinprev, tmp1, &u, &udagger);
			}
			else if (mu == 2){
				dslash_2(spinnext, spinprev, tmp1, &u, &udagger);
			}
			else{
				dslash_3(spinnext, spinprev, tmp1, &u, &udagger);
			}
		
			//second two
			get_spinor_from_field(in, spinnext, nextplusnu_spat, nextplusnu_time);
			get_spinor_from_field(in, spinprev, nextminusnu_spat, nextminusnu_time);

		//TODO
// 		if(prev_time == 0) spinor_apply_bc(spinnext, theta);
// 		else if(prevminusnu_time == NTIME-1) spinor_apply_bc(spinprev, theta);
      
			get_su3matrix(&u,gaugefield,next_spat, next_time,nu[i]);
			get_su3matrix(&udagger,gaugefield,nextminusnu_spat,nextminusnu_time,mu);
			adjoin_su3matrix(&udagger);

			//update links with chemical potential, this shall be put into compiler option lateron
			gaugefield_apply_chem_pot(&u, &udagger, chem_pot_re, chem_pot_im);

			//store in tmp2
			if(mu == 0){
				dslash_0(spinnext, spinprev, tmp2, &u, &udagger);
			}
			else if (mu == 1){
				dslash_1(spinnext, spinprev, tmp2, &u, &udagger);
			}
			else if (mu == 2){
				dslash_2(spinnext, spinprev, tmp2, &u, &udagger);
			}
			else{
				dslash_3(spinnext, spinprev, tmp2, &u, &udagger);
			}
		}
	
		get_su3matrix(&u,gaugefield,prev_spat, prev_time,mu);
		get_su3matrix(&udagger,gaugefield,pos,t,mu);
		adjoin_su3matrix(&udagger);
		
		if(mu == 0){
// 			set_local_zero_spinor(tmp1);
// 			//add spinors to get the sum nu neq mu
// 			spinors_accumulate(tmp1, tmp4);
// 			spinors_accumulate(tmp1, tmp5);
// 			spinors_accumulate(tmp1, tmp6);
// 			spinors_accumulate(tmp2, tmp8);
// 			spinors_accumulate(tmp2, tmp9);
// 			spinors_accumulate(tmp2, tmp10);
			ddaggerslash_0(tmp1, tmp2, spinout, &u, &udagger);
		}
		else if (mu == 1){
			ddaggerslash_1(tmp1, tmp2, spinout, &u, &udagger);
		}
		else if (mu == 2){
			ddaggerslash_2(tmp1, tmp2, spinout, &u, &udagger);
		}
		else{
			ddaggerslash_3(tmp1, tmp2, spinout, &u, &udagger);
		}
	}
	return;
}

