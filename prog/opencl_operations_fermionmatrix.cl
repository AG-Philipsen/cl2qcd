/** @file
 * Device code for operations on the fermion matrix
 */

//opencl_fermionmatrix.cl

#ifdef _FERMIONS_
void inline dslash_spatial (hmc_spinor * spinout, int * coord, int dir, int pos, int t, __global hmc_spinor_field* in,  __global hmc_ocl_gaugefield* gaugefield, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im){

	int next, prev;
	hmc_spinor spinnext[SPINORSIZE];
	hmc_spinor spinprev[SPINORSIZE];
	hmc_ocl_su3matrix u[SU3SIZE];
	hmc_ocl_su3matrix udagger[SU3SIZE]; 

	next = get_neighbor(pos,dir);
	prev = get_lower_neighbor(pos,dir);
  
	get_spinor_from_field(in, spinnext, next, t);
	get_spinor_from_field(in, spinprev, prev, t);
      
	get_su3matrix(u,gaugefield,pos,t,dir);
	get_su3matrix(udagger,gaugefield,prev,t,dir);
	adjoin_su3matrix(udagger);

#ifdef _NPBC_S_
	if(coord[dir] == NSPACE-1) spinor_apply_bc(spinnext, theta);
	else if(coord[dir] == 0) spinor_apply_bc(spinprev, theta);
#endif

	if(dir == 1) dslash_1(spinnext, spinprev, spinout, u, udagger);
	else if(dir == 2) dslash_2(spinnext, spinprev, spinout, u, udagger);
	else dslash_3(spinnext, spinprev, spinout, u, udagger);

	return;
}

void inline dslash_temporal (hmc_spinor * spinout, int pos, int t, __global hmc_spinor_field* in,  __global hmc_ocl_gaugefield* gaugefield, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im){
	int next, prev;
	hmc_spinor spinnext[SPINORSIZE];
	hmc_spinor spinprev[SPINORSIZE];
	hmc_ocl_su3matrix u[SU3SIZE];
	hmc_ocl_su3matrix udagger[SU3SIZE]; 

	next = (t+1)%NTIME; 
	prev = (t-1+NTIME)%NTIME;

	get_spinor_from_field(in, spinnext, pos, next);
	get_spinor_from_field(in, spinprev, pos, prev);

#ifdef _NPBC_T_
	if(next == 0) spinor_apply_bc(spinnext, theta);
	else if(prev == NTIME-1) spinor_apply_bc(spinprev, theta);
#endif

	get_su3matrix(u,gaugefield,pos,t,0);
	get_su3matrix(udagger,gaugefield,pos,prev,0);
	adjoin_su3matrix(udagger);

#ifdef _CP_REAL_
	gaugefield_apply_chem_pot_real(u, udagger, chem_pot_re);
#endif
#ifdef _CP_IMAG_
	gaugefield_apply_chem_pot_imag(u, udagger, chem_pot_im);
#endif
      
	dslash_0(spinnext, spinprev, spinout, u, udagger);

	return;
}

void inline dslash_spatial_eoprec (hmc_spinor * spinout, int * coord, int dir, int pos, int t, __global hmc_eoprec_spinor_field* in,  __global hmc_ocl_gaugefield* gaugefield, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im){

	int next, prev;
	int neo_next;
	int neo_prev;
	hmc_spinor spinnext[SPINORSIZE];
	hmc_spinor spinprev[SPINORSIZE];
	hmc_ocl_su3matrix u[SU3SIZE];
	hmc_ocl_su3matrix udagger[SU3SIZE]; 

	next = get_neighbor(pos,dir);
	prev = get_lower_neighbor(pos,dir);
  
	neo_next = get_n_eoprec(next, t);
	neo_prev = get_n_eoprec(prev, t);

	get_spinor_from_eoprec_field(in,spinnext,neo_next);
	get_spinor_from_eoprec_field(in,spinprev,neo_prev);
	
	get_su3matrix(u,gaugefield,pos,t,dir);
	get_su3matrix(udagger,gaugefield,prev,t,dir);
	adjoin_su3matrix(udagger);

#ifdef _NPBC_S_
	if(coord[dir] == NSPACE-1) spinor_apply_bc(spinnext, theta);
	else if(coord[dir] == 0) spinor_apply_bc(spinprev, theta);
#endif

	if(dir == 1) dslash_1(spinnext, spinprev, spinout, u, udagger);
	else if(dir == 2) dslash_2(spinnext, spinprev, spinout, u, udagger);
	else dslash_3(spinnext, spinprev, spinout, u, udagger);

	return;
}

void inline dslash_temporal_eoprec (hmc_spinor * spinout, int pos, int t, __global hmc_eoprec_spinor_field* in,  __global hmc_ocl_gaugefield* gaugefield, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im){
	int next, prev;
	int neo_next;
	int neo_prev;
	hmc_spinor spinnext[SPINORSIZE];
	hmc_spinor spinprev[SPINORSIZE];
	hmc_ocl_su3matrix u[SU3SIZE];
	hmc_ocl_su3matrix udagger[SU3SIZE]; 

	next = (t+1)%NTIME; 
	prev = (t-1+NTIME)%NTIME;

	neo_next = get_n_eoprec(pos, next);
	neo_prev = get_n_eoprec(pos, prev);

	get_spinor_from_eoprec_field(in,spinnext,neo_next);
	get_spinor_from_eoprec_field(in,spinprev,neo_prev);

#ifdef _NPBC_T_
	if(next == 0) spinor_apply_bc(spinnext, theta);
	else if(prev == NTIME-1) spinor_apply_bc(spinprev, theta);
#endif

	get_su3matrix(u,gaugefield,pos,t,0);
	get_su3matrix(udagger,gaugefield,pos,prev,0);
	adjoin_su3matrix(udagger);

#ifdef _CP_REAL_
	gaugefield_apply_chem_pot_real(u, udagger, chem_pot_re);
#endif
#ifdef _CP_IMAG_
	gaugefield_apply_chem_pot_imag(u, udagger, chem_pot_im);
#endif
      
	dslash_0(spinnext, spinprev, spinout, u, udagger);

	return;
}

__kernel void M_diag (__global hmc_spinor_field* in, __global hmc_spinor_field* out) {
	int id = get_global_id(0);
	int id_tmp;
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id(0);

	int pos;
	int t;
	hmc_float kappa_tmp = KAPPA;
	hmc_float mu_tmp = MU;
	hmc_spinor spinout[SPINORSIZE];
	for(id_tmp = id; id_tmp < VOL4D/2; id_tmp += global_size){
		get_even_site(id_tmp, &pos, &t);
		get_spinor_from_field(in,spinout,pos,t);
		M_diag_local(spinout, kappa_tmp, mu_tmp);
		put_spinor_to_field(spinout,out,pos,t);

		get_odd_site(id_tmp, &pos, &t);
		get_spinor_from_field(in,spinout,pos,t);
		M_diag_local(spinout, kappa_tmp, mu_tmp);
		put_spinor_to_field(spinout,out,pos,t);
	}

	return;
}

__kernel void M_sitediagonal (__global hmc_eoprec_spinor_field* in, __global hmc_eoprec_spinor_field* out) {
	int id = get_global_id(0);
	int id_tmp;
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id(0);

	hmc_float kappa_tmp = KAPPA;
	hmc_float mu_tmp = MU;
	hmc_spinor spinout[SPINORSIZE];
	for(id_tmp = id; id_tmp < VOL4D/2; id_tmp += global_size){
		get_spinor_from_eoprec_field(in,spinout,id_tmp);
		M_diag_local(spinout, kappa_tmp, mu_tmp);
		put_spinor_to_eoprec_field(spinout,out,id_tmp);
	}

	return;
}

__kernel void M_inverse_sitediagonal (__global hmc_eoprec_spinor_field* in, __global hmc_eoprec_spinor_field* out) {
	int id = get_global_id(0);
	int id_tmp;
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id(0);

	hmc_float kappa_tmp = KAPPA;
	hmc_float mu_tmp = MU;
	hmc_spinor spinout[SPINORSIZE];
	for(id_tmp = id; id_tmp < VOL4D/2; id_tmp += global_size){
		hmc_float minuskappa = -kappa_tmp;
		get_spinor_from_eoprec_field(in,spinout,id_tmp);
		M_diag_local(spinout, minuskappa, mu_tmp);
		hmc_float denom = 1. + 4.*kappa_tmp*kappa_tmp*mu_tmp*mu_tmp;
		real_multiply_spinor(spinout,1./denom);
		put_spinor_to_eoprec_field(spinout,out,id_tmp);
	}

	return;
}

__kernel void dslash(__global hmc_spinor_field* in, __global hmc_spinor_field* out, __global hmc_ocl_gaugefield* gaugefield, __global hmc_float* theta_in, __global hmc_float* chem_pot_re_in, __global hmc_float* chem_pot_im_in){
	int id = get_global_id(0);
	int global_size = get_global_size(0);
	int id_tmp;
	int local_size = get_local_size(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id(0);

	int pos;
	int t;

	hmc_float theta = *theta_in;
	hmc_float chem_pot_im = *chem_pot_im_in;
	hmc_float chem_pot_re = *chem_pot_re_in;
	
	hmc_spinor spinout[SPINORSIZE];

	for(id_tmp = id; id_tmp < VOL4D/2; id_tmp += global_size){
		
		//CP: this can become NDIM-1 for mem-saving
		int coord[NDIM];
		coord[0]=0;		

		set_local_zero_spinor(spinout);
		get_even_site(id_tmp, &pos, &t);
		for(int j=1;j<NDIM;j++) coord[j] = get_spacecoord(pos,j);
		
		// spinout = U_0*(r-gamma_0)*spinnext + U^dagger_0(x-hat0) * (r+gamma_0)*spinprev
		dslash_temporal (spinout, pos, t, in, gaugefield, theta, chem_pot_re, chem_pot_im);
		// spinout += U_1*(r-gamma_1)*spinnext + U^dagger_1(x-hat1) * (r+gamma_1)*spinprev 
		dslash_spatial (spinout, coord, 1, pos, t, in, gaugefield, theta, chem_pot_re, chem_pot_im);
		// spinout += U_2*(r-gamma_2)*spinnext + U^dagger_2(x-hat2) * (r+gamma_2)*spinprev
		dslash_spatial (spinout, coord, 2, pos, t, in, gaugefield, theta, chem_pot_re, chem_pot_im);
		// spinout += U_3*(r-gamma_3)*spinnext + U^dagger_3(x-hat3) * (r+gamma_3)*spinprev
		dslash_spatial (spinout, coord, 3, pos, t, in, gaugefield, theta, chem_pot_re, chem_pot_im);
  
		put_spinor_to_field(spinout,out,pos,t);

		set_local_zero_spinor(spinout);
		get_odd_site(id_tmp, &pos, &t);
		for(int j=1;j<NDIM;j++) coord[j] = get_spacecoord(pos,j);
		
		// spinout = U_0*(r-gamma_0)*spinnext + U^dagger_0(x-hat0) * (r+gamma_0)*spinprev
		dslash_temporal (spinout, pos, t, in, gaugefield, theta, chem_pot_re, chem_pot_im);
		// spinout += U_1*(r-gamma_1)*spinnext + U^dagger_1(x-hat1) * (r+gamma_1)*spinprev 
		dslash_spatial (spinout, coord, 1, pos, t, in, gaugefield, theta, chem_pot_re, chem_pot_im);
		// spinout += U_2*(r-gamma_2)*spinnext + U^dagger_2(x-hat2) * (r+gamma_2)*spinprev
		dslash_spatial (spinout, coord, 2, pos, t, in, gaugefield, theta, chem_pot_re, chem_pot_im);
		// spinout += U_3*(r-gamma_3)*spinnext + U^dagger_3(x-hat3) * (r+gamma_3)*spinprev
		dslash_spatial (spinout, coord, 3, pos, t, in, gaugefield, theta, chem_pot_re, chem_pot_im);
    
		put_spinor_to_field(spinout,out,pos,t);
	}
	return;
}

__kernel void dslash_eoprec(__global hmc_eoprec_spinor_field* in, __global hmc_eoprec_spinor_field* out, __global hmc_ocl_gaugefield* gaugefield, __global hmc_float* theta_in, __global hmc_float* chem_pot_re_in, __global hmc_float* chem_pot_im_in, int evenodd){
	int id = get_global_id(0);
	int global_size = get_global_size(0);
	int id_tmp;
	int local_size = get_local_size(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id(0);

	int pos;
	int t;

	hmc_float theta = *theta_in;
	hmc_float chem_pot_im = *chem_pot_im_in;
	hmc_float chem_pot_re = *chem_pot_re_in;
	hmc_float kappa = KAPPA;

	hmc_spinor spinout[SPINORSIZE];

	for(id_tmp = id; id_tmp < VOL4D/2; id_tmp += global_size){
		
		//CP: this can become NDIM-1 for mem-saving
		int coord[NDIM];
		coord[0]=0;		

		if(evenodd == ODD) get_odd_site(id_tmp, &pos, &t);
    else get_even_site(id_tmp, &pos, &t);
		
		set_local_zero_spinor(spinout);
		for(int j=1;j<NDIM;j++) coord[j] = get_spacecoord(pos,j);
		
		// spinout = U_0*(r-gamma_0)*spinnext + U^dagger_0(x-hat0) * (r+gamma_0)*spinprev
		dslash_temporal_eoprec (spinout, pos, t, in, gaugefield, theta, chem_pot_re, chem_pot_im);
		// spinout += U_1*(r-gamma_1)*spinnext + U^dagger_1(x-hat1) * (r+gamma_1)*spinprev 
		dslash_spatial_eoprec (spinout, coord, 1, pos, t, in, gaugefield, theta, chem_pot_re, chem_pot_im);
		// spinout += U_2*(r-gamma_2)*spinnext + U^dagger_2(x-hat2) * (r+gamma_2)*spinprev
		dslash_spatial_eoprec (spinout, coord, 2, pos, t, in, gaugefield, theta, chem_pot_re, chem_pot_im);
		// spinout += U_3*(r-gamma_3)*spinnext + U^dagger_3(x-hat3) * (r+gamma_3)*spinprev
		dslash_spatial_eoprec (spinout, coord, 3, pos, t, in, gaugefield, theta, chem_pot_re, chem_pot_im);
  
		real_multiply_spinor(spinout,-kappa);
		put_spinor_to_eoprec_field(spinout,out,id_tmp);
	}
	return;
}

#endif //_FERMIONS_
