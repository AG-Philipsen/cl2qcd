/**
 @file fermionmatrix-functions for eoprec spinorfields
*/

//"local" dslash working on a particular link (n,t) of an eoprec field
//NOTE: each component is multiplied by +KAPPA, so the resulting spinor has to be mutliplied by -1 to obtain the correct dslash!!!
//the difference to the "normal" dslash is that the coordinates of the neighbors have to be transformed into an eoprec index
spinor inline dslash_eoprec_local_0(__global const hmc_float * const restrict in, __global hmc_float const * const restrict field, const int n, const int t)
{
	spinor out_tmp, plus;
	int dir, nn, nn_eo;
	su3vec psi, phi;
	Matrixsu3 U;
	//this is used to save the BC-conditions...
	hmc_complex bc_tmp;
	out_tmp = set_spinor_zero();

	//go through the different directions
	///////////////////////////////////
	// mu = 0
	///////////////////////////////////
	dir = 0;
	///////////////////////////////////
	//mu = +0
	nn = (t + 1) % NTIME;
	//transform normal indices to eoprec index
	nn_eo = get_n_eoprec(n, nn);
	plus = getSpinorSOA_eo(in, nn_eo);
	U = getSU3SOA(field, get_global_link_pos_SOA(dir, n, t));
	//if chemical potential is activated, U has to be multiplied by appropiate factor
#ifdef _CP_REAL_
	U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
	hmc_complex cpi_tmp = {COSCPI, SINCPI};
	U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
	bc_tmp.re = KAPPA_TEMPORAL_RE;
	bc_tmp.im = KAPPA_TEMPORAL_IM;
	///////////////////////////////////
	// Calculate psi/phi = (1 - gamma_0) plus/y
	// with 1 - gamma_0:
	// | 1  0  1  0 |        | psi.e0 + psi.e2 |
	// | 0  1  0  1 |  psi = | psi.e1 + psi.e3 |
	// | 1  0  1  0 |        | psi.e1 + psi.e3 |
	// | 0  1  0  1 |        | psi.e0 + psi.e2 |
	///////////////////////////////////
	// psi = 0. component of (1-gamma_0)y
	psi = su3vec_acc(plus.e0, plus.e2);
	// phi = U*psi
	phi =  su3matrix_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
	out_tmp.e2 = su3vec_acc(out_tmp.e2, psi);
	// psi = 1. component of (1-gamma_0)y
	psi = su3vec_acc(plus.e1, plus.e3);
	// phi = U*psi
	phi =  su3matrix_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
	out_tmp.e3 = su3vec_acc(out_tmp.e3, psi);

	/////////////////////////////////////
	//mu = -0
	nn = (t - 1 + NTIME) % NTIME;
	//transform normal indices to eoprec index
	nn_eo = get_n_eoprec(n, nn);
	plus = getSpinorSOA_eo(in, nn_eo);
	U = getSU3SOA(field, get_global_link_pos_SOA(dir, n, nn));
	//if chemical potential is activated, U has to be multiplied by appropiate factor
	//this is the same as at mu=0 in the imag. case, since U is taken to be U^+ later:
	//  (exp(iq)U)^+ = exp(-iq)U^+
	//as it should be
	//in the real case, one has to take exp(q) -> exp(-q)
#ifdef _CP_REAL_
	U = multiply_matrixsu3_by_real (U, MEXPCPR);
#endif
#ifdef _CP_IMAG_
	hmc_complex cpi_tmp = {COSCPI, SINCPI};
	U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
	//in direction -mu, one has to take the complex-conjugated value of bc_tmp. this is done right here.
	bc_tmp.re = KAPPA_TEMPORAL_RE;
	bc_tmp.im = MKAPPA_TEMPORAL_IM;
	///////////////////////////////////
	// Calculate psi/phi = (1 + gamma_0) y
	// with 1 + gamma_0:
	// | 1  0 -1  0 |       | psi.e0 - psi.e2 |
	// | 0  1  0 -1 | psi = | psi.e1 - psi.e3 |
	// |-1  0  1  0 |       | psi.e1 - psi.e2 |
	// | 0 -1  0  1 |       | psi.e0 - psi.e3 |
	///////////////////////////////////
	// psi = 0. component of (1+gamma_0)y
	psi = su3vec_dim(plus.e0, plus.e2);
	// phi = U*psi
	phi = su3matrix_dagger_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e0 = su3vec_acc(out_tmp.e0, psi);
	out_tmp.e2 = su3vec_dim(out_tmp.e2, psi);
	// psi = 1. component of (1+gamma_0)y
	psi = su3vec_dim(plus.e1, plus.e3);
	// phi = U*psi
	phi = su3matrix_dagger_times_su3vec(U, psi);
	psi = su3vec_times_complex(phi, bc_tmp);
	out_tmp.e1 = su3vec_acc(out_tmp.e1, psi);
	out_tmp.e3 = su3vec_dim(out_tmp.e3, psi);

	return out_tmp;
}

spinor inline dslash_eoprec_local_1(__global const hmc_float * const restrict in, __global hmc_float  const * const restrict field, const int n, const int t)
{
	spinor out_tmp, plus;
	int dir, nn, nn_eo;
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
	//transform normal indices to eoprec index
	nn_eo = get_n_eoprec(nn, t);
	plus = getSpinorSOA_eo(in, nn_eo);
	U = getSU3SOA(field, get_global_link_pos_SOA(dir, n, t));
	bc_tmp.re = KAPPA_SPATIAL_RE;
	bc_tmp.im = KAPPA_SPATIAL_IM;
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
	//transform normal indices to eoprec index
	nn_eo = get_n_eoprec(nn, t);
	plus = getSpinorSOA_eo(in, nn_eo);
	U = getSU3SOA(field, get_global_link_pos_SOA(dir, nn, t));
	//in direction -mu, one has to take the complex-conjugated value of bc_tmp. this is done right here.
	bc_tmp.re = KAPPA_SPATIAL_RE;
	bc_tmp.im = MKAPPA_SPATIAL_IM;
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

spinor inline dslash_eoprec_local_2(__global const hmc_float * const restrict in, __global hmc_float const * const restrict field, const int n, const int t)
{
	spinor out_tmp, plus;
	int dir, nn, nn_eo;
	su3vec psi, phi;
	Matrixsu3 U;
	//this is used to save the BC-conditions...
	hmc_complex bc_tmp;
	out_tmp = set_spinor_zero();

	///////////////////////////////////
	// mu = 2
	///////////////////////////////////
	dir = 2;

	///////////////////////////////////
	// mu = +2
	nn = get_neighbor(n, dir);
	//transform normal indices to eoprec index
	nn_eo = get_n_eoprec(nn, t);
	plus = getSpinorSOA_eo(in, nn_eo);
	U = getSU3SOA(field, get_global_link_pos_SOA(dir, n, t));
	bc_tmp.re = KAPPA_SPATIAL_RE;
	bc_tmp.im = KAPPA_SPATIAL_IM;
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
	//transform normal indices to eoprec index
	nn_eo = get_n_eoprec(nn, t);
	plus = getSpinorSOA_eo(in, nn_eo);
	U = getSU3SOA(field, get_global_link_pos_SOA(dir, nn, t));
	//in direction -mu, one has to take the complex-conjugated value of bc_tmp. this is done right here.
	bc_tmp.re = KAPPA_SPATIAL_RE;
	bc_tmp.im = MKAPPA_SPATIAL_IM;
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

spinor inline dslash_eoprec_local_3(__global const hmc_float * const restrict in, __global hmc_float const * const restrict field, const int n, const int t)
{
	spinor out_tmp, plus;
	int dir, nn, nn_eo;
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
	//transform normal indices to eoprec index
	nn_eo = get_n_eoprec(nn, t);
	plus = getSpinorSOA_eo(in, nn_eo);
	U = getSU3SOA(field, get_global_link_pos_SOA(dir, n, t));
	bc_tmp.re = KAPPA_SPATIAL_RE;
	bc_tmp.im = KAPPA_SPATIAL_IM;
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
	//transform normal indices to eoprec index
	nn_eo = get_n_eoprec(nn, t);
	plus = getSpinorSOA_eo(in, nn_eo);
	U = getSU3SOA(field, get_global_link_pos_SOA(dir, nn, t));
	//in direction -mu, one has to take the complex-conjugated value of bc_tmp. this is done right here.
	bc_tmp.re = KAPPA_SPATIAL_RE;
	bc_tmp.im = MKAPPA_SPATIAL_IM;
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

