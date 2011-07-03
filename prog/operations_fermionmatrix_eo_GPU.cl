/**
 @file fermionmatrix-functions for eoprec spinorfields
*/

//"local" dslash working on a particular link (n,t) of an eoprec field
//NOTE: each component is multiplied by +KAPPA, so the resulting spinor has to be mutliplied by -1 to obtain the correct dslash!!!
//the difference to the "normal" dslash is that the coordinates of the neighbors have to be transformed into an eoprec index
spinor dslash_eoprec_local_0(__global spinorfield_eoprec * in,__global ocl_s_gaugefield * field, int n, int t){
	spinor out_tmp, plus;
	int dir, nn, nn_eo;
	su3vec psi, phi;
	Matrixsu3 U;
	hmc_float kappa_minus = KAPPA;
	
	out_tmp = set_spinor_zero();

	//go through the different directions
	///////////////////////////////////
	// mu = 0
	///////////////////////////////////
	dir = 0;
	///////////////////////////////////
	//mu = +0
	nn = (t+1)%NTIME;
	//transform normal indices to eoprec index
	nn_eo = get_n_eoprec(n, nn);
	plus = get_spinor_from_eoprec_field(in, nn_eo);
	U = field[get_global_link_pos(dir, n, t)];
	///////////////////////////////////
	// Calculate psi/phi = (1 - gamma_0) y
	// with 1 - gamma_0:
	// | 1  0 -1  0 |        |       psi.e0 - psi.e2 |
	// | 0  1  0 -1 |  psi = |       psi.e1 - psi.e3 |
	// |-1  0  1  0 |        |(-1)*( psi.e1 - psi.e3 |
	// | 0 -1  0  1 |        |(-1)*( psi.e0 - psi.e2 |
	///////////////////////////////////
	// psi = 0. component of (1-gamma_0)y
	psi = su3vec_dim(plus.e0, plus.e2);
	// phi = U*psi
	phi =  su3matrix_times_su3vec(U, psi);
	psi = su3vec_times_real(phi, kappa_minus);
	out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
	out_tmp.e2 = su3vec_dim(out_tmp.e2, psi);
	// psi = 1. component of (1-gamma_0)y
	psi = su3vec_dim(plus.e1, plus.e3);
	// phi = U*psi
	phi =  su3matrix_times_su3vec(U, psi);
	psi = su3vec_times_real(phi, kappa_minus);
	out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
	out_tmp.e3 = su3vec_dim(out_tmp.e3, psi);

	/////////////////////////////////////
	//mu = -0
	nn = (t-1+NTIME)%NTIME;
	//transform normal indices to eoprec index
	nn_eo = get_n_eoprec(n, nn);
	plus = get_spinor_from_eoprec_field(in, nn_eo);
	U = field[get_global_link_pos(dir, n, nn)];
	///////////////////////////////////
	// Calculate psi/phi = (1 + gamma_0) y
	// with 1 + gamma_0:
	// | 1  0  1  0 |       | psi.e0 + i*psi.e2 |
	// | 0  1  0  1 | psi = | psi.e1 + i*psi.e3 |
	// | 1  0  1  0 |       | psi.e1 + i*psi.e2 |
	// | 0  1  0  1 |       | psi.e0 + i*psi.e3 |
	///////////////////////////////////
	// psi = 0. component of (1+gamma_0)y
	psi = su3vec_acc(plus.e0, plus.e2);
	// phi = U*psi
	phi = su3matrix_dagger_times_su3vec(U, psi);
	psi = su3vec_times_real(phi, kappa_minus);
	out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
	out_tmp.e2 = su3vec_acc(out_tmp.e2, psi);
	// psi = 1. component of (1+gamma_0)y
	psi = su3vec_acc(plus.e1, plus.e3);
	// phi = U*psi
	phi = su3matrix_dagger_times_su3vec(U, psi);
	psi = su3vec_times_real(phi, kappa_minus);
	out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
	out_tmp.e3 = su3vec_acc(out_tmp.e3, psi);		

	return out_tmp;
}
	
	
spinor dslash_eoprec_local_1(__global spinorfield_eoprec * in,__global ocl_s_gaugefield * field, int n, int t){
	spinor out_tmp, plus;
	int dir, nn, nn_eo;
	su3vec psi, phi;
	Matrixsu3 U;
	hmc_float kappa_minus = KAPPA;
	
	out_tmp = set_spinor_zero();
	//CP: all actions correspond to the mu = 0 ones
	///////////////////////////////////
	// mu = 1
	///////////////////////////////////
	dir = 1;
	
	///////////////////////////////////
	// mu = +1
	nn = get_neighbor(n,dir);
	//transform normal indices to eoprec index
	nn_eo = get_n_eoprec(nn, t);
	plus = get_spinor_from_eoprec_field(in, nn_eo);
	U = field[get_global_link_pos(dir, n, t)];
	///////////////////////////////////
	// Calculate (1 - gamma_1) y
	// with 1 - gamma_1:
	// | 1  0  0 -i |       |      psi.e0 - i*psi.e3  |
	// | 0  1 -i  0 | psi = |      psi.e1 - i*psi.e2  |
	// | 0  i  1  0 |       |(i)*( psi.e1 - i*psi.e2) |
	// | i  0  0  1 |       |(i)*( psi.e0 - i*psi.e3) |
	///////////////////////////////////
	psi = su3vec_dim_i(plus.e0, plus.e3);
	phi = su3matrix_times_su3vec(U, psi);
	psi = su3vec_times_real(phi, kappa_minus);
	out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
	out_tmp.e3 = su3vec_acc_i(out_tmp.e3, psi);
	
	psi = su3vec_dim_i(plus.e1, plus.e2);
	phi = su3matrix_times_su3vec(U, psi);
	psi = su3vec_times_real(phi, kappa_minus);
	out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
	out_tmp.e2 = su3vec_acc_i(out_tmp.e2, psi);		

	///////////////////////////////////
	//mu = -1
	nn = get_lower_neighbor(n,dir);
	//transform normal indices to eoprec index
	nn_eo = get_n_eoprec(nn, t);
	plus = get_spinor_from_eoprec_field(in, nn_eo);
	U = field[get_global_link_pos(dir, nn, t)];
	///////////////////////////////////
	// Calculate (1 + gamma_1) y
	// with 1 + gamma_1:
	// | 1  0  0  i |       |       psi.e0 + i*psi.e3  |
	// | 0  1  i  0 | psi = |       psi.e1 + i*psi.e2  |
	// | 0 -i  1  0 |       |(-i)*( psi.e1 + i*psi.e2) |
	// |-i  0  0  1 |       |(-i)*( psi.e0 + i*psi.e3) |
	///////////////////////////////////
	psi = su3vec_acc_i(plus.e0, plus.e3);
	phi = su3matrix_dagger_times_su3vec(U, psi);
	psi = su3vec_times_real(phi, kappa_minus);
	out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
	out_tmp.e3 = su3vec_dim_i(out_tmp.e3, psi);
	
	psi = su3vec_acc_i(plus.e1, plus.e2);
	phi = su3matrix_dagger_times_su3vec(U, psi);
	psi = su3vec_times_real(phi, kappa_minus);
	out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
	out_tmp.e2 = su3vec_dim_i(out_tmp.e2, psi);		

	return out_tmp;
}
	
	
spinor dslash_eoprec_local_2(__global spinorfield_eoprec * in,__global ocl_s_gaugefield * field, int n, int t){
	spinor out_tmp, plus;
	int dir, nn, nn_eo;
	su3vec psi, phi;
	Matrixsu3 U;
	hmc_float kappa_minus = KAPPA;
	
	out_tmp = set_spinor_zero();
	///////////////////////////////////
	// mu = 2
	///////////////////////////////////
	dir = 2;
	
	///////////////////////////////////
	// mu = +2
	nn = get_neighbor(n,dir);
	//transform normal indices to eoprec index
	nn_eo = get_n_eoprec(nn, t);
	plus = get_spinor_from_eoprec_field(in, nn_eo);
	U = field[get_global_link_pos(dir, n, t)];
	///////////////////////////////////
	// Calculate (1 - gamma_2) y
	// with 1 - gamma_2:
	// | 1  0  0  1 |       |       psi.e0 + psi.e3  |
	// | 0  1 -1  0 | psi = |       psi.e1 - psi.e2  |
	// | 0 -1  1  0 |       |(-1)*( psi.e1 + psi.e2) |
	// | 1  0  0  1 |       |     ( psi.e0 + psi.e3) |
	///////////////////////////////////
	psi = su3vec_acc(plus.e0, plus.e3);
	phi = su3matrix_times_su3vec(U, psi);
	psi = su3vec_times_real(phi, kappa_minus);
	out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
	out_tmp.e3 = su3vec_acc(out_tmp.e3, psi);
	
	psi = su3vec_dim(plus.e1, plus.e2);
	phi = su3matrix_times_su3vec(U, psi);
	psi = su3vec_times_real(phi, kappa_minus);
	out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
	out_tmp.e2 = su3vec_dim(out_tmp.e2, psi);	
	
	///////////////////////////////////
	//mu = -2
	nn = get_lower_neighbor(n,dir);
	//transform normal indices to eoprec index
	nn_eo = get_n_eoprec(nn, t);
	plus = get_spinor_from_eoprec_field(in, nn_eo);
	U = field[get_global_link_pos(dir, nn, t)];
	///////////////////////////////////
	// Calculate (1 + gamma_2) y
	// with 1 + gamma_2:
	// | 1  0  0 -1 |       |       psi.e0 - psi.e3  |
	// | 0  1  1  0 | psi = |       psi.e1 + psi.e2  |
	// | 0  1  1  0 |       |     ( psi.e1 + psi.e2) |
	// |-1  0  0  1 |       |(-1)*( psi.e0 - psi.e3) |
	///////////////////////////////////
	psi = su3vec_dim(plus.e0, plus.e3);
	phi = su3matrix_dagger_times_su3vec(U, psi);
	psi = su3vec_times_real(phi, kappa_minus);
	out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
	out_tmp.e3 = su3vec_dim(out_tmp.e3, psi);
	
	psi = su3vec_acc(plus.e1, plus.e2);
	phi = su3matrix_dagger_times_su3vec(U, psi);
	psi = su3vec_times_real(phi, kappa_minus);
	out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
	out_tmp.e2 = su3vec_acc(out_tmp.e2, psi);	
	return out_tmp;
}

	
spinor dslash_eoprec_local_3(__global spinorfield_eoprec * in,__global ocl_s_gaugefield * field, int n, int t){
	spinor out_tmp, plus;
	int dir, nn, nn_eo;
	su3vec psi, phi;
	Matrixsu3 U;
	hmc_float kappa_minus = KAPPA;
	
	out_tmp = set_spinor_zero();
	///////////////////////////////////
	// mu = 3
	///////////////////////////////////
	dir = 3;
	
	///////////////////////////////////
	// mu = +3
	nn = get_neighbor(n,dir);
	//transform normal indices to eoprec index
	nn_eo = get_n_eoprec(nn, t);
	plus = get_spinor_from_eoprec_field(in, nn_eo);
	U = field[get_global_link_pos(dir, n, t)];
	///////////////////////////////////
	// Calculate (1 - gamma_3) y
	// with 1 - gamma_3:
	// | 1  0 -i  0 |        |       psi.e0 - i*psi.e2  |
	// | 0  1  0  i |  psi = |       psi.e1 + i*psi.e3  |
	// | i  0  1  0 |        |   i *(psi.e0 - i*psi.e2) |
	// | 0 -i  0  1 |        | (-i)*(psi.e1 + i*psi.e3) |
	///////////////////////////////////
	psi = su3vec_dim_i(plus.e0, plus.e2);
	phi = su3matrix_times_su3vec(U, psi);
	psi = su3vec_times_real(phi, kappa_minus);
	out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
	out_tmp.e2 = su3vec_acc_i(out_tmp.e2, psi);
	
	psi = su3vec_acc_i(plus.e1, plus.e3);
	phi = su3matrix_times_su3vec(U, psi);
	psi = su3vec_times_real(phi, kappa_minus);
	out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
	out_tmp.e3 = su3vec_dim_i(out_tmp.e3, psi);	

	///////////////////////////////////
	//mu = -3
	nn = get_lower_neighbor(n,dir);
	//transform normal indices to eoprec index
	nn_eo = get_n_eoprec(nn, t);
	plus = get_spinor_from_eoprec_field(in, nn_eo);
	U = field[get_global_link_pos(dir, nn, t)];
	///////////////////////////////////
	// Calculate (1 + gamma_3) y
	// with 1 + gamma_3:
	// | 1  0  i  0 |       |       psi.e0 + i*psi.e2  |
	// | 0  1  0 -i | psi = |       psi.e1 - i*psi.e3  |
	// |-i  0  1  0 |       | (-i)*(psi.e0 + i*psi.e2) |
	// | 0  i  0  1 |       |   i *(psi.e1 - i*psi.e3) |
	///////////////////////////////////
	psi = su3vec_acc_i(plus.e0, plus.e2);
	phi = su3matrix_dagger_times_su3vec(U, psi);
	psi = su3vec_times_real(phi, kappa_minus);
	out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
	out_tmp.e2 = su3vec_dim_i(out_tmp.e2, psi);
	
	psi = su3vec_dim_i(plus.e1, plus.e3);
	phi = su3matrix_dagger_times_su3vec(U, psi);
	psi = su3vec_times_real(phi, kappa_minus);
	out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
	out_tmp.e3 = su3vec_acc_i(out_tmp.e3, psi);

	return out_tmp;
}

__kernel void dslash_eoprec(__global spinorfield_eoprec* in, __global spinorfield_eoprec* out, __global ocl_s_gaugefield* field, int evenodd){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	int n,t;
	spinor out_tmp;
	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE; id_tmp += global_size) {	
		if(evenodd == ODD) get_odd_site(id_tmp, &n, &t);
		else get_even_site(id_tmp, &n, &t);
	
		out_tmp = dslash_eoprec_local_0(in, field, n, t);
		out_tmp = spinor_acc(out_tmp, dslash_eoprec_local_1(in, field, n, t));
		out_tmp = spinor_acc_acc(out_tmp, dslash_eoprec_local_2(in, field, n, t),dslash_eoprec_local_3(in, field, n, t));
		put_spinor_to_eoprec_field(out_tmp, out, id_tmp);

	}
}

__kernel void M_sitediagonal(__global spinorfield_eoprec * in, __global spinorfield_eoprec * out){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	spinor out_tmp;
	spinor plus;
	
	hmc_complex twistfactor = {1., MUBAR};
	hmc_complex twistfactor_minus = {1., MMUBAR};
	
	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE; id_tmp += global_size) {	
		out_tmp = set_spinor_zero();
		//get input spinor
		plus = get_spinor_from_eoprec_field(in, id_tmp);
		//Diagonalpart:
		out_tmp = M_diag_local(plus, twistfactor, twistfactor_minus);
		put_spinor_to_eoprec_field(out_tmp, out, id_tmp);
	}
}

__kernel void M_inverse_sitediagonal(__global spinorfield_eoprec * in, __global spinorfield_eoprec * out){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	spinor out_tmp;
	spinor plus;
	hmc_complex twistfactor = {1., MUBAR};
	hmc_complex twistfactor_minus = {1., MMUBAR};
	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE; id_tmp += global_size) {	
		out_tmp = set_spinor_zero();
		//get input spinor
		plus = get_spinor_from_eoprec_field(in, id_tmp);
		//Diagonalpart, here the twisted factor give the inverse matrix:
		out_tmp = M_diag_local(plus, twistfactor_minus, twistfactor);
		hmc_float denom = 1./(1. + MUBAR*MUBAR);
		out_tmp = real_multiply_spinor(out_tmp,denom);
		put_spinor_to_eoprec_field(out_tmp, out, id_tmp);
	}
}

__kernel void gamma5_eoprec(__global spinorfield_eoprec *in, __global spinorfield_eoprec *out){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	spinor out_tmp;

	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE; id_tmp += global_size) {	

		out_tmp = get_spinor_from_eoprec_field(in, id_tmp); 
		out_tmp = gamma5_local(out_tmp);
		put_spinor_to_eoprec_field(out_tmp, out, id_tmp);
	}
}
