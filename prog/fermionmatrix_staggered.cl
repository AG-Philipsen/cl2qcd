/**
 @file staggered fermionmatrix-functions that are called from within kernels
*/

/*
 * These should go into an own file...
 */

su3vec get_su3vec_from_field(__global const su3vec * const restrict in, const int n, const int t)
{
	int pos = get_global_pos(n, t);
	su3vec out;
	out = in[pos];
	return out;
}

void put_su3vec_to_field(const su3vec in, __global su3vec * const restrict out, const int n, const int t)
{
	int pos = get_global_pos(n, t);
	out[pos] = in;
}

//"local" dslash working on a particular link (n,t) in a specific direction
//spinor dslash_local_0(__global const spinorfield * const restrict in,__global const ocl_s_gaugefield * const restrict field, const int n, const int t){
su3vec dslash_local_0(__global const su3vec * const restrict in, __global const Matrixsu3StorageType * const restrict field, int n, int t)
{
	su3vec out_tmp, plus;
	int dir, nn;
	Matrixsu3 U;
	//this is used to save the BC-conditions...
	hmc_complex bc_tmp;
	out_tmp = set_su3vec_zero();
	//this is used to save the staggered phase...
	hmc_float st_phase = 1.;

	//go through the different directions
	///////////////////////////////////
	// mu = 0
	///////////////////////////////////
	dir = 0;
	///////////////////////////////////
	//mu = +0
	nn = get_neighbor_temporal(t);
	plus = get_su3vec_from_field(in, n, nn);
	U = getSU3(field, get_link_pos(dir, n, t));
	//if chemical potential is activated, U has to be multiplied by appropiate factor
#ifdef _CP_REAL_
	U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
	hmc_complex cpi_tmp = {COSCPI, SINCPI};
	U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
	bc_tmp.re = TEMPORAL_RE;
	bc_tmp.im = TEMPORAL_IM;
	///////////////////////////////////
	// chi = U*plus
	out_tmp =  su3matrix_times_su3vec(U, plus);
	out_tmp = su3vec_times_complex(out_tmp, bc_tmp);

	/////////////////////////////////////
	//mu = -0
	nn = get_lower_neighbor_temporal(t);
	plus = get_su3vec_from_field(in, n, nn);
	U = getSU3(field, get_link_pos(dir, n, nn));
	//if chemical potential is activated, U has to be multiplied by appropiate factor
	//this is the same as at mu=0 in the imag. case, since U is taken to be U^+ later:
	//  (exp(iq)U)^+ = exp(-iq)U^+
	//as it should be
	//in the real case, one has to take exp(q) -> exp(-q)
#ifdef _CP_REAL_
	U = multiply_matrixsu3_by_real (U, MEXPCPR);
#endif
#ifdef _CP_IMAG_
	hmc_complex cpi_tmp2 = {COSCPI, SINCPI};
	U = multiply_matrixsu3_by_complex (U, cpi_tmp2 );
#endif
	//in direction -mu, one has to take the complex-conjugated value of bc_tmp. this is done right here.
	bc_tmp.re = TEMPORAL_RE;
	bc_tmp.im = MTEMPORAL_IM;
	// chi = U^dagger*plus
	out_tmp = su3matrix_dagger_times_su3vec(U, plus);
	out_tmp = su3vec_times_complex(out_tmp, bc_tmp);

	//multiply with staggered phase
	out_tmp = su3vec_times_real(out_tmp, st_phase);

	return out_tmp;
}

//these have to be adjusted according to the above function..
/*
spinor dslash_local_1(__global const spinor * const restrict in, __global const Matrixsu3StorageType * const restrict field, int n, int t, hmc_float kappa_in)
{
	spinor out_tmp, plus;
	int dir, nn;
	su3vec psi, phi;
	Matrixsu3 U;
	//this is used to save the BC-conditions...
	hmc_complex bc_tmp;
	out_tmp = set_spinor_zero();

	//CP: all actions correspond to the mu = 0 ones
	///////////////////////////////////
	// mu = 1
	///////////////////////////////////
	dir = 1;

	///////////////////////////////////
	// mu = +1
	nn = get_neighbor(n, dir);
	plus = get_spinor_from_field(in, nn, t);
	U = getSU3(field, get_link_pos(dir, n, t));
	bc_tmp.re = kappa_in * SPATIAL_RE;
	bc_tmp.im = kappa_in * SPATIAL_IM;
	/////////////////////////////////
	//Calculate (1 - gamma_1) y
	//with 1 - gamma_1:
	//| 1  0  0  i |       |       psi.e0 + i*psi.e3  |
	//| 0  1  i  0 | psi = |       psi.e1 + i*psi.e2  |
	//| 0  i  1  0 |       |(-i)*( psi.e1 + i*psi.e2) |
	//| i  0  0  1 |       |(-i)*( psi.e0 + i*psi.e3) |
	/////////////////////////////////
	psi = su3vec_acc_i(plus.e0, plus.e3);
	phi = su3matrix_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
	out_tmp.e3 = su3vec_dim_i(out_tmp.e3, psi);

	psi = su3vec_acc_i(plus.e1, plus.e2);
	phi = su3matrix_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
	out_tmp.e2 = su3vec_dim_i(out_tmp.e2, psi);

	///////////////////////////////////
	//mu = -1
	nn = get_lower_neighbor(n, dir);
	plus = get_spinor_from_field(in, nn, t);
	U = getSU3(field, get_link_pos(dir, nn, t));
	//in direction -mu, one has to take the complex-conjugated value of bc_tmp. this is done right here.
	bc_tmp.re = kappa_in * SPATIAL_RE;
	bc_tmp.im = kappa_in * MSPATIAL_IM;
	///////////////////////////////////
	// Calculate (1 + gamma_1) y
	// with 1 + gamma_1:
	// | 1  0  0 -i |       |       psi.e0 - i*psi.e3  |
	// | 0  1 -i  0 | psi = |       psi.e1 - i*psi.e2  |
	// | 0  i  1  0 |       |(-i)*( psi.e1 - i*psi.e2) |
	// | i  0  0  1 |       |(-i)*( psi.e0 - i*psi.e3) |
	///////////////////////////////////
	psi = su3vec_dim_i(plus.e0, plus.e3);
	phi = su3matrix_dagger_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
	out_tmp.e3 = su3vec_acc_i(out_tmp.e3, psi);

	psi = su3vec_dim_i(plus.e1, plus.e2);
	phi = su3matrix_dagger_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
	out_tmp.e2 = su3vec_acc_i(out_tmp.e2, psi);


	return out_tmp;
}

spinor dslash_local_2(__global const spinor * const restrict in, __global const Matrixsu3StorageType * const restrict field, int n, int t, hmc_float kappa_in)
{
	spinor out_tmp, plus;
	int dir, nn;
	su3vec psi, phi;
	Matrixsu3 U;
	//this is used to save the BC-conditions...
	hmc_complex bc_tmp;
	out_tmp = set_spinor_zero();;

	///////////////////////////////////
	// mu = 2
	///////////////////////////////////
	dir = 2;

	///////////////////////////////////
	// mu = +2
	nn = get_neighbor(n, dir);
	plus = get_spinor_from_field(in, nn, t);
	U = getSU3(field, get_link_pos(dir, n, t));
	bc_tmp.re = kappa_in * SPATIAL_RE;
	bc_tmp.im = kappa_in * SPATIAL_IM;
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
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
	out_tmp.e3 = su3vec_acc(out_tmp.e3, psi);

	psi = su3vec_dim(plus.e1, plus.e2);
	phi = su3matrix_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
	out_tmp.e2 = su3vec_dim(out_tmp.e2, psi);

	///////////////////////////////////
	//mu = -2
	nn = get_lower_neighbor(n, dir);
	plus = get_spinor_from_field(in,  nn, t);
	U = getSU3(field, get_link_pos(dir, nn, t));
	//in direction -mu, one has to take the complex-conjugated value of bc_tmp. this is done right here.
	bc_tmp.re = kappa_in * SPATIAL_RE;
	bc_tmp.im = kappa_in * MSPATIAL_IM;
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
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
	out_tmp.e3 = su3vec_dim(out_tmp.e3, psi);

	psi = su3vec_acc(plus.e1, plus.e2);
	phi = su3matrix_dagger_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
	out_tmp.e2 = su3vec_acc(out_tmp.e2, psi);

	return out_tmp;
}

spinor dslash_local_3(__global const spinor * const restrict in, __global const Matrixsu3StorageType * const restrict field, int n, int t, hmc_float kappa_in)
{
	spinor out_tmp, plus;
	int dir, nn;
	su3vec psi, phi;
	Matrixsu3 U;
	//this is used to save the BC-conditions...
	hmc_complex bc_tmp;
	out_tmp = set_spinor_zero();

	///////////////////////////////////
	// mu = 3
	///////////////////////////////////
	dir = 3;

	///////////////////////////////////
	// mu = +3
	nn = get_neighbor(n, dir);
	plus = get_spinor_from_field(in, nn, t);
	U = getSU3(field, get_link_pos(dir, n, t));
	bc_tmp.re = kappa_in * SPATIAL_RE;
	bc_tmp.im = kappa_in * SPATIAL_IM;
	///////////////////////////////////
	// Calculate (1 - gamma_3) y
	// with 1 - gamma_3:
	// | 1  0  i  0 |        |       psi.e0 + i*psi.e2  |
	// | 0  1  0 -i |  psi = |       psi.e1 - i*psi.e3  |
	// |-i  0  1  0 |        |   i *(psi.e0 + i*psi.e2) |
	// | 0  i  0  1 |        | (-i)*(psi.e1 - i*psi.e3) |
	///////////////////////////////////
	psi = su3vec_acc_i(plus.e0, plus.e2);
	phi = su3matrix_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
	out_tmp.e2 = su3vec_dim_i(out_tmp.e2, psi);

	psi = su3vec_dim_i(plus.e1, plus.e3);
	phi = su3matrix_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
	out_tmp.e3 = su3vec_acc_i(out_tmp.e3, psi);

	///////////////////////////////////
	//mu = -3
	nn = get_lower_neighbor(n, dir);
	plus = get_spinor_from_field(in, nn, t);
	U = getSU3(field, get_link_pos(dir, nn, t));
	//in direction -mu, one has to take the complex-conjugated value of bc_tmp. this is done right here.
	bc_tmp.re = kappa_in * SPATIAL_RE;
	bc_tmp.im = kappa_in * MSPATIAL_IM;
	///////////////////////////////////
	// Calculate (1 + gamma_3) y
	// with 1 + gamma_3:
	// | 1  0 -i  0 |       |       psi.e0 - i*psi.e2  |
	// | 0  1  0  i | psi = |       psi.e1 + i*psi.e3  |
	// | i  0  1  0 |       | (-i)*(psi.e0 - i*psi.e2) |
	// | 0 -i  0  1 |       |   i *(psi.e1 + i*psi.e3) |
	///////////////////////////////////
	psi = su3vec_dim_i(plus.e0, plus.e2);
	phi = su3matrix_dagger_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
	out_tmp.e2 = su3vec_acc_i(out_tmp.e2, psi);

	psi = su3vec_acc_i(plus.e1, plus.e3);
	phi = su3matrix_dagger_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
	out_tmp.e3 = su3vec_dim_i(out_tmp.e3, psi);

	return out_tmp;
}
*/
