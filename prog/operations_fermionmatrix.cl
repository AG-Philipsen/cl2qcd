/**
 @file fermionmatrix
*/

// void print_su3vec(su3vec in){
//      printf("(%f,%f)\t(%f,%f)\t(%f,%f)\t", in.e0.re, in.e0.im, in.e1.re, in.e1.im, in.e2.re, in.e2.im);
// }
// 
// void print_spinor(spinor in){
//      print_su3vec(in.e0);
//      print_su3vec(in.e1);
//      print_su3vec(in.e2);
//      print_su3vec(in.e3);
//      printf("\n");
// }


__kernel void M(__global spinorfield * in, __global ocl_s_gaugefield * field, __global spinorfield * out){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	int n,t,nn, dir;
	spinor out_tmp;
	spinor plus;
	Matrixsu3 U;
	su3vec psi, phi;
	
	hmc_complex twistfactor = {1., MUBAR};
	hmc_complex twistfactor_minus = {1., MMUBAR};
	hmc_float kappa_minus = MKAPPA;

	/** @todo implement BC! Beware, the kappa at -mu then has to be complex conjugated (see tmlqcd)*/
	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size) {	
	
		/** @todo this must be done more efficient */
		if(id_tmp%2 == 0) get_even_site(id_tmp/2, &n, &t);
		else get_odd_site(id_tmp/2, &n, &t);
		
		out_tmp = set_spinor_zero();
		//get input spinor
		plus = get_spinor_from_field(in, n, t);
		//Diagonalpart:
		//	(1+i*mubar*gamma_5)psi = (1, mubar)psi.0,1 (1,-mubar)psi.2,3
		out_tmp.e0 = su3vec_times_complex(plus.e0, twistfactor);
		out_tmp.e1 = su3vec_times_complex(plus.e1, twistfactor);
		out_tmp.e2 = su3vec_times_complex(plus.e2, twistfactor_minus);
		out_tmp.e3 = su3vec_times_complex(plus.e3, twistfactor_minus);
		
		//go through the different directions
		///////////////////////////////////
		// mu = 0
		///////////////////////////////////
		//dir = 0;
		dir = 1;
		///////////////////////////////////
		//mu = +0
		nn = (t+1)%NTIME;
		plus = get_spinor_from_field(in, n, nn);
		U = unit_matrixsu3();	
		//U = field[get_global_link_pos(dir, n, t)];
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
		plus = get_spinor_from_field(in, n, nn);
		U = unit_matrixsu3();
		//U = field[get_global_link_pos(dir, n, nn)];
		//get_su3matrix(&U, field, n, nn, dir);
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

		//CP: all actions correspond to the mu = 0 ones
		///////////////////////////////////
		// mu = 1
		///////////////////////////////////
		dir = 1;
		
		///////////////////////////////////
		// mu = +1
		nn = get_neighbor(n,dir);
		plus = get_spinor_from_field(in, nn, t);
		U = unit_matrixsu3();
		//U = field[get_global_link_pos(dir, n, t)];
		//get_su3matrix(&U, field, n, t, dir);
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
		plus = get_spinor_from_field(in, nn, t);
		U = unit_matrixsu3();
		//U = field[get_global_link_pos(dir, nn, t)];
		//get_su3matrix(&U, field, nn, t, dir);
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
		
		///////////////////////////////////
		// mu = 2
		///////////////////////////////////
		dir = 2;
		
		///////////////////////////////////
		// mu = +2
		nn = get_neighbor(n,dir);
		plus = get_spinor_from_field(in, nn, t);
		U = unit_matrixsu3();
		//U = field[get_global_link_pos(dir, n, t)];
		//get_su3matrix(&U, field, n, t, dir);
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
		plus = get_spinor_from_field(in,  nn, t);
		U = unit_matrixsu3();
		//U = field[get_global_link_pos(dir, nn, t)];
		//get_su3matrix(&U, field, nn, t, dir);
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

		///////////////////////////////////
		// mu = 3
		///////////////////////////////////
		dir = 3;
		
		///////////////////////////////////
		// mu = +3
		nn = get_neighbor(n,dir);
		plus = get_spinor_from_field(in, nn, t);
		U = unit_matrixsu3();
		//U = field[get_global_link_pos(dir, n, t)];
		//get_su3matrix(&U, field, n, t, dir);
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
		plus = get_spinor_from_field(in, nn, t);
		U = unit_matrixsu3();
		//U = field[get_global_link_pos(dir, nn, t)];
		//get_su3matrix(&U, field, nn, t, dir);
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
	
		put_spinor_to_field(out_tmp, out, n, t);
	}
}

//CP: perhaps these are still needed
__kernel void M_diag(){

}

__kernel void dslash(){


}

__kernel void dslash_eoprec(){


}

__kernel void M_sitediagonal(){

}

__kernel void M_inverse_sitediagonal(){

}
