

/*
void inline dslash_spatial (hmc_spinor * spinout, int * coord, int dir, int pos, int t, __global hmc_spinor_field* in,  __global hmc_gaugefield* gaugefield, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im){

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
	gaugefield_apply_chem_pot(u, udagger, chem_pot_re, chem_pot_im);

	if(coord[dir] == NSPACE-1) spinor_apply_bc(spinnext, theta);
	else if(coord[dir] == 0) spinor_apply_bc(spinprev, theta);
      
	multiply_spinor_gamma1(spinnext,tmp);
	real_multiply_spinor(tmp,-hmc_one_f);
	spinors_accumulate(spinnext,tmp);
	su3matrix_times_spinor(&u,spinnext,tmp);
	spinors_accumulate(spinout,tmp);

	multiply_spinor_gamma1(spinprev,tmp);
	spinors_accumulate(spinprev,tmp);
	su3matrix_times_spinor(&u,spinprev,tmp);
	spinors_accumulate(spinout,tmp);

	return;
}

void inline dslash_temporal (hmc_spinor * spinout, int * coord, int pos, int t, __global hmc_spinor_field* in,  __global hmc_gaugefield* gaugefield, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im){
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
      
	get_su3matrix(u,gaugefield,pos,t,0);
	get_su3matrix(udagger,gaugefield,pos,prev,0);
	adjoin_su3matrix(udagger);

	//update links with chemical potential, this shall be put into compiler option lateron
	gaugefield_apply_chem_pot(u, udagger, chem_pot_re, chem_pot_im);
      
	multiply_spinor_gamma0(spinnext,tmp);
	real_multiply_spinor(tmp,-hmc_one_f);
	spinors_accumulate(spinnext,tmp);
	su3matrix_times_spinor(&u,spinnext,tmp);
	spinors_accumulate(spinout,tmp);

	multiply_spinor_gamma0(spinprev,tmp);
	spinors_accumulate(spinprev,tmp);
	su3matrix_times_spinor(&udagger,spinprev,tmp);
	spinors_accumulate(spinout,tmp);

	return;
}

//!! the for loop
//!! for(id_tmp = id; id_tmp < VOL4D/2/num_groups*(group_id + 1 ); id_tmp += local_size)
//!! can also be
//!! for(id_tmp=id; id_tmp< VOL4D/2; id_tmp+=get_num_groups(0)*local_size)
//!! for the beginning, I stick to the first version!!

__kernel void M_diag (__global hmc_spinor_field* in, __global hmc_spinor_field* out, hmc_float kappa, hmc_float mu) {
	int id = get_global_id(0);
	int id_tmp;
	int local_size = get_local_size(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id(0);

	int pos;
	int t;
	hmc_spinor spinout[SPINORSIZE];
	hmc_spinor spintmp[SPINORSIZE];

	for(id_tmp = id; id_tmp < VOL4D/2/num_groups*(group_id + 1 ); id_tmp += local_size){
		get_even_side(id_tmp, &pos, &t);
		get_spinor_from_field(in,spinout,pos,t);
		hmc_float twistfactor = 2*kappa*mu;
		multiply_spinor_factor_gamma5(spinout,spintmp,twistfactor);
		spinors_accumulate(spinout,spintmp);
		put_spinor_to_field(spinout,out,pos,t);

		get_odd_side(id_tmp, &pos, &t);
		get_spinor_from_field(in,spinout,pos,t);
		hmc_float twistfactor = 2*kappa*mu;
		multiply_spinor_factor_gamma5(spinout,spintmp,twistfactor);
		spinors_accumulate(spinout,spintmp);
		put_spinor_to_field(spinout,out,pos,t);
	}

	return;
}

__kernel void dslash(__global hmc_spinor_field* in, __global hmc_spinor_field* out, __global hmc_gaugefield* gaugefield, hmc_float theta, hmc_float chem_pot_re, hmc_float chem_pot_im){
	int id = get_global_id(0);
	int id_tmp;
	int local_size = get_local_size(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id(0);

	int pos;
	int t;

	hmc_spinor spinout[SPINORSIZE];
	hmc_spinor spinnext[SPINORSIZE];
	hmc_spinor spinprev[SPINORSIZE];
	hmc_spinor tmp[SPINORSIZE];
	hmc_su3matrix u;
	hmc_su3matrix udagger;
	int next;
	int prev;

	//CP: this is done without odd sites only for space-saving
	for(int id_tmp = id; id_tmp < VOL4D/num_groups*(group_id + 1 ); id_tmp += local_size){
		set_local_zero_spinor(spinout);
		get_even_side (id_tmp, &pos, &t);
		//CP: this can become NDIM-1 for mem-saving
		int coord[NDIM];
		coord[0]=0;
		for(int j=1;j<NDIM;j++) coord[j] = get_spacecoord(pos,j);
		
		// spinout = U_0*(r-gamma_0)*spinnext + U^dagger_0(x-hat0) * (r+gamma_0)*spinprev
		dslash_temporal (spinout, coord, pos, t, in, gaugefield, theta, chem_pot_re, chem_pot_im);
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

// -alpha*x + y
//CP: defined with a minus!!!
__kernel void saxpy(__global hmc_spinor_field* x, __global hmc_spinor_field* y, hmc_float * alpha, __global hmc_spinor_field* out){
	int id = get_global_id(0);
	int id_tmp;
	int local_size = get_local_size(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id(0);

	int pos;
	int t;
	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE/num_groups*(group_id + 1 ); id_tmp += local_size){
		hmc_complex tmp1 = complexmult(alpha,&x[id_tmp]);
		((out)[id_tmp]).re = -(tmp1).re + y[id_tmp].re;
		((out)[id_tmp]).im = -(tmp1).im + y[id_tmp].im;
	}

	return;
}

//alpha*x + beta*y + z
__kernel void saxsbypz(__global hmc_spinor_field* x, __global hmc_spinor_field* y, __global hmc_spinor_field* z, hmc_float * alpha, hmc_float * beta, __global hmc_spinor_field* out){
	int id = get_global_id(0);
	int id_tmp;
	int local_size = get_local_size(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id(0);

	int pos;
	int t;
	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE/num_groups*(group_id + 1 ); id_tmp += local_size){
		hmc_complex tmp1 = complexmult(alpha,&x[id_tmp]);
    hmc_complex tmp2 = complexmult(beta,&y[id_tmp]);
    ((out)[id_tmp]).re = (tmp1).re + (tmp2).re + z[id_tmp].re;
    ((out)[id_tmp]).im = (tmp1).im + (tmp2).im + z[id_tmp].im;
	}

	return;
}

// complex (!!!) scalarproduct, return in result 
// --> use 2 kernels: 1 for the summation in one block and 1 for summation over blockresults
__kernel void scalar_product( __global __readonly hmc_spinor_field *x, __global __readonly hmc_spinor_field *y, __global __writeonly hmc_complex* result, __local hmc_complex* result_local ){
	int local_size = get_local_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	hmc_complex sum;
	sum.re = 0.;
	sum.im = 0.;

	//!! CP: perhaps here one can first copy a whole spinor and then do the dot-prod
	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE/num_groups*(group_id + 1 ); id_tmp += local_size){
    sum.re += x[id_tmp].re*y[id_tmp].re + x[id_tmp].im*y[id_tmp].im;
    sum.im += x[id_tmp].re*y[id_tmp].im - x[id_tmp].im*y[id_tmp].re;
	}
	
	// sync threads
	barrier();
	//reduction
	if (id<64)
		(result_local[id]).re=sum.re;
		(result_local[id]).im=sum.im;
	barrier();
	if (id>=64)
		complex_accumulate(sum, result_local[id-64]);
//!! CP: this must be result_local[id%64]+=sum; ??
	barrier();
	if (id<32)
		(result_local[id]).re=sum.re;
		(result_local[id]).im=sum.im;
	barrier();
	if (threadid.x>=32)
		complex_accumulate(sum, result_local[id-32]);
	barrier();
	if (id<16)
		(result_local[id]).re=sum.re;
		(result_local[id]).im=sum.im;
	barrier();
	if (threadid.x>=16)
		complex_accumulate(sum, result_local[id-16]);
	barrier();
	if (id<8)
		(result_local[id]).re=sum.re;
		(result_local[id]).im=sum.im;
	barrier();
	//thread 0 sums up the result_local and stores it in array result
	if (id==0)
		for (int i = 0; i<8; i++) complex_accumulate(result[get_group_id(0)], result_local[i]);

	return;
}

__kernel void scalar_product_reduction(__global __readonly hmc_float* result_tmp, __global __writeonly hmc_float* result){
	for (int i = 1; i<get_num_groups(0); i++){
		complex_accumulate(result_tmp[0], result_tmp[i]);
	}
	//!!CP: pointer????
	(*result) = result_tmp[0];

	return;
}

__kernel void copy_spinorfield(const hmc_spinorfield *x, const hmc_spinorfield *y){
	int local_size = get_local_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);

	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE/num_groups*(group_id + 1 ); id_tmp += local_size){
		y[id_tmp] = x[id_tmp];
	}
	return;
}

__kernel void set_zero_spinorfield(const hmc_spinorfield *x){
	int local_size = get_local_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);

	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE/num_groups*(group_id + 1 ); id_tmp += local_size){
		x[id_tmp] = 0.;
	}
	return;
}

// complex (!!!) squarenorm, return in result 
// --> use 2 kernels: 1 for the summation in one block and 1 for summation over blockresults (this is the same as for the scalar_product)
__kernel void global_squarenorm( __global __readonly hmc_spinor_field *x, __global __writeonly hmc_float* result, __local hmc_float* result_local ){
	int local_size = get_local_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	hmc_float sum;
	sum = 0.;

	//!! CP: perhaps here one can first copy a whole spinor and then do the dot-prod
	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE/num_groups*(group_id + 1 ); id_tmp += local_size){
    sum += x[id_tmp].re*x[id_tmp].re + x[id_tmp].im*x[id_tmp].im;
	}
	
	// sync threads
	barrier();
	//reduction
	if (id<64)
		(result_local[id])=sum;
	barrier();
	if (id>=64)
		result_local[id-64]+=sum;
//!! CP: this must be result_local[id%64]+=sum; ??
	barrier();
	if (id<32)
		(result_local[id])=sum;
	barrier();
	if (threadid.x>=32)
		result_local[id-32]+=sum;
	barrier();
	if (id<16)
		(result_local[id])=sum;
	barrier();
	if (threadid.x>=16)
		result_local[id-16]+=sum;
	barrier();
	if (id<8)
		(result_local[id])=sum;
	barrier();
	//thread 0 sums up the result_local and stores it in array result
	if (id==0)
		result[ get_group_id(0) ] = 	result_local[0] + result_local[1] + result_local[2] + result_local[3] + 
				                result_local[4] + result_local[5] + result_local[6] + result_local[7];
	return;
}

__kernel void scalar_product_reduction(__global __readonly hmc_float* result_tmp, __global __writeonly hmc_float* result){
	for (int i = 1; i<get_num_groups(0); i++){
		complex_accumulate(result_tmp[0], result_tmp[i]);
	}
	//!!CP: pointer????
	(*result) = result_tmp[0];

	return;
}

__kernel void ratio(__global __readonly hmc_complex * a, __global __readonly hmc_complex * b, __global __writeonly hmc_complex * out){
	(*out) =  complexdivide(a, b);
	return;
}

__kernel void product(__global __readonly hmc_complex * a, __global __readonly hmc_complex * b, __global __writeonly hmc_complex * out){
	(*out) =  complexmult(a, b);
	return;
}

__kernel void set_zero_spinorfield( __global __writeonly hmc_spinor_field *x ){
	int local_size = get_local_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);

	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE/num_groups*(group_id + 1 ); id_tmp += local_size){
		x[id_tmp].re = 0.;
		x[id_tmp].im = 0.;
	}

*/