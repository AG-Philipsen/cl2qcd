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

//local Diagonalmatrix:
//	(1+i*mubar*gamma_5)psi = (1, mubar)psi.0,1 (1,-mubar)psi.2,3
spinor inline M_diag_local(spinor in, hmc_complex factor1, hmc_complex factor2){
	spinor tmp;
	tmp.e0 = su3vec_times_complex(in.e0, factor1);
	tmp.e1 = su3vec_times_complex(in.e1, factor1);
	tmp.e2 = su3vec_times_complex(in.e2, factor2);
	tmp.e3 = su3vec_times_complex(in.e3, factor2);
	return tmp;	
}

/** @todo this can be optimized... */
//local gamma5:
//	(gamma_5)psi = (1)psi.0,1 (-1)psi.2,3
spinor inline gamma5_local(spinor in){
	spinor tmp;
	tmp.e0 = in.e0;
	tmp.e1 = in.e1;
	tmp.e2 = su3vec_times_real(in.e2, -1.);
	tmp.e3 = su3vec_times_real(in.e3, -1.);
	return tmp;	
}

//"local" dslash working on a particular link (n,t)
//NOTE: each component is multiplied by +KAPPA, so the resulting spinor has to be mutliplied by -1 to obtain the correct dslash!!!
spinor inline dslash_local(__global spinorfield * in,__global ocl_s_gaugefield * field, int n, int t){
	spinor out_tmp, plus;
	int dir, nn;
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
	plus = get_spinor_from_field(in, n, nn);
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
	plus = get_spinor_from_field(in, n, nn);
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

	//CP: all actions correspond to the mu = 0 ones
	///////////////////////////////////
	// mu = 1
	///////////////////////////////////
	dir = 1;
	
	///////////////////////////////////
	// mu = +1
	nn = get_neighbor(n,dir);
	plus = get_spinor_from_field(in, nn, t);
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
	plus = get_spinor_from_field(in, nn, t);
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
	
	///////////////////////////////////
	// mu = 2
	///////////////////////////////////
	dir = 2;
	
	///////////////////////////////////
	// mu = +2
	nn = get_neighbor(n,dir);
	plus = get_spinor_from_field(in, nn, t);
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
	plus = get_spinor_from_field(in,  nn, t);
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

	///////////////////////////////////
	// mu = 3
	///////////////////////////////////
	dir = 3;
	
	///////////////////////////////////
	// mu = +3
	nn = get_neighbor(n,dir);
	plus = get_spinor_from_field(in, nn, t);
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
	plus = get_spinor_from_field(in, nn, t);
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

//"local" dslash working on a particular link (n,t) of an eoprec field
//NOTE: each component is multiplied by +KAPPA, so the resulting spinor has to be mutliplied by -1 to obtain the correct dslash!!!
//the difference to the "normal" dslash is that the coordinates of the neighbors have to be transformed into an eoprec index
spinor inline dslash_eoprec_local(__global spinorfield_eoprec * in,__global ocl_s_gaugefield * field, int n, int t){
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

__kernel void M(__global spinorfield * in, __global ocl_s_gaugefield * field, __global spinorfield * out){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	int n,t;
	spinor out_tmp, out_tmp2;
	spinor plus;
	
	hmc_complex twistfactor = {1., MUBAR};
	hmc_complex twistfactor_minus = {1., MMUBAR};

	/** @todo implement BC! Beware, the kappa at -mu then has to be complex conjugated (see tmlqcd)*/
	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size) {	
	
		/** @todo this must be done more efficient */
		if(id_tmp%2 == 0) get_even_site(id_tmp/2, &n, &t);
		else get_odd_site(id_tmp/2, &n, &t);
		
		out_tmp = set_spinor_zero();
		//get input spinor
		plus = get_spinor_from_field(in, n, t);
		//Diagonalpart:
		out_tmp = M_diag_local(plus, twistfactor, twistfactor_minus);
		//calc dslash (this includes mutliplication with kappa)
		out_tmp2 = dslash_local(in, field, n, t);
		//M = M_diag - dslash
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		put_spinor_to_field(out_tmp, out, n, t);
	}
}

//CP: perhaps these are still needed
__kernel void M_diag(__global spinorfield * in, __global spinorfield * out){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	int n,t;
	spinor out_tmp;
	spinor plus;
	
	hmc_complex twistfactor = {1., MUBAR};
	hmc_complex twistfactor_minus = {1., MMUBAR};
	
	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size) {	
		
		/** @todo this must be done more efficient */
		if(id_tmp%2 == 0) get_even_site(id_tmp/2, &n, &t);
		else get_odd_site(id_tmp/2, &n, &t);
		
		out_tmp = set_spinor_zero();
		//get input spinor
		plus = get_spinor_from_field(in, n, t);
		//Diagonalpart:
		out_tmp = M_diag_local(plus, twistfactor, twistfactor_minus);
		put_spinor_to_field(out_tmp, out, n, t);
	}
}

//to be removed....
__kernel void dslash(){
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
		out_tmp = dslash_eoprec_local(in, field, n, t);
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

__kernel void gamma5(__global spinorfield *in, __global spinorfield *out){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	spinor out_tmp;
	int n,t;

	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size) {	
		/** @todo this must be done more efficient */
		if(id_tmp%2 == 0) get_even_site(id_tmp/2, &n, &t);
		else get_odd_site(id_tmp/2, &n, &t);

		out_tmp = get_spinor_from_field(in,n, t); 
		out_tmp = gamma5_local(out_tmp);
		put_spinor_to_field(out_tmp, out, n, t);
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

__kernel void Qplus(__global spinorfield * in, __global ocl_s_gaugefield * field, __global spinorfield * out){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	int n,t;
	spinor out_tmp, out_tmp2;
	spinor plus;
	
	hmc_complex twistfactor = {1., MUBAR};
	hmc_complex twistfactor_minus = {1., MMUBAR};

	/** @todo implement BC! Beware, the kappa at -mu then has to be complex conjugated (see tmlqcd)*/
	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size) {	
	
		/** @todo this must be done more efficient */
		if(id_tmp%2 == 0) get_even_site(id_tmp/2, &n, &t);
		else get_odd_site(id_tmp/2, &n, &t);
		
		out_tmp = set_spinor_zero();
		//get input spinor
		plus = get_spinor_from_field(in, n, t);
		//Diagonalpart: this is just the normal tm-diagonal matrix
		out_tmp = M_diag_local(plus, twistfactor, twistfactor_minus);
		//calc dslash (this includes mutliplication with kappa)
		out_tmp2 = dslash_local(in, field, n, t);
		//M = M_diag - dslash
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		//gamma_5
		out_tmp = gamma5_local(out_tmp);
		put_spinor_to_field(out_tmp, out, n, t);
	}
}

__kernel void Qminus(__global spinorfield * in, __global ocl_s_gaugefield * field, __global spinorfield * out){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	int n,t;
	spinor out_tmp, out_tmp2;
	spinor plus;
	
	hmc_complex twistfactor = {1., MUBAR};
	hmc_complex twistfactor_minus = {1., MMUBAR};

	/** @todo implement BC! Beware, the kappa at -mu then has to be complex conjugated (see tmlqcd)*/
	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size) {	
	
		/** @todo this must be done more efficient */
		if(id_tmp%2 == 0) get_even_site(id_tmp/2, &n, &t);
		else get_odd_site(id_tmp/2, &n, &t);
		
		out_tmp = set_spinor_zero();
		//get input spinor
		plus = get_spinor_from_field(in, n, t);
		//Diagonalpart: this is normal tm-diagonal matrix with negative imaginary part
		out_tmp = M_diag_local(plus, twistfactor_minus, twistfactor);
		//calc dslash (this includes mutliplication with kappa)
		out_tmp2 = dslash_local(in, field, n, t);
		//M = M_diag - dslash
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		//gamma_5
		out_tmp = gamma5_local(out_tmp);
		put_spinor_to_field(out_tmp, out, n, t);
	}
}