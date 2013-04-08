/**
 * @file kernel for the eoprec fermion force
 */

//Here, two cases are possible:
//evenodd = ODD or EVEN
//	EVEN corresponds to case where an odd site of Y is connected to an even site of X (because if x is odd x+mu is even)
//	ODD is just the other way around.
//The difference to the non-eo kernel is that one has to calculate the eo position of the "plus" spinor out of the neighbour coordinates on each occasion. This is done just like in the dslash_eo kernel!
__kernel void fermion_force_eo(__global const Matrixsu3StorageType * const restrict field, __global const spinorStorageType * const restrict Y, __global const spinorStorageType * const restrict X, __global aeStorageType * const restrict out, int evenodd, hmc_float kappa_in)
{
	// must include HALO, as we are updating neighbouring sites
	// -> not all local sites will fully updated if we don't calculate on halo indices, too
	PARALLEL_FOR(id_mem, EOPREC_SPINORFIELDSIZE_MEM) {
		//caculate (pos,time) out of id_local depending on evenodd
		//st_index pos = (evenodd == ODD) ? get_even_st_idx_local(id_local) : get_odd_st_idx_local(id_local);
		st_index pos = (evenodd == ODD) ? get_even_st_idx(id_mem) : get_odd_st_idx(id_mem);

		Matrixsu3 U;
		Matrix3x3 v1, v2, tmp;
		su3vec psia, psib, phia, phib;
		spinor y, plus;
		int nn;
		ae out_tmp;
		int global_link_pos;
		int global_link_pos_down;
		//this is used to save the BC-conditions...
		hmc_complex bc_tmp;
		int dir;
		int n = pos.space;
		int t = pos.time;
		int nn_eo;

		y = getSpinor_eo(Y, get_eo_site_idx_from_st_idx(pos));
		///////////////////////////////////
		// Calculate gamma_5 y
		///////////////////////////////////
		y = gamma5_local(y);

		//go through the different directions
		///////////////////////////////////
		// mu = 0
		///////////////////////////////////
		dir = 0;
		//the 2 here comes from Tr(lambda_ij) = 2delta_ij
		bc_tmp.re = 2.* kappa_in * TEMPORAL_RE;
		bc_tmp.im = 2.* kappa_in * TEMPORAL_IM;

		///////////////////////////////////
		//mu = +0
		global_link_pos = get_link_pos(dir, n, t);
		nn = get_neighbor_temporal(t);
		//transform normal indices to eoprec index
		nn_eo = get_n_eoprec(n, nn);
		plus = getSpinor_eo(X, nn_eo);
		U = get_matrixsu3(field, n, t, dir);
		//if chemical potential is activated, U has to be multiplied by appropiate factor
#ifdef _CP_REAL_
		U = multiply_matrixsu3_by_real (U, EXPCPR);
#endif
#ifdef _CP_IMAG_
		hmc_complex cpi_tmp = {COSCPI, SINCPI};
		U = multiply_matrixsu3_by_complex (U, cpi_tmp );
#endif
		///////////////////////////////////
		// Calculate psi/phi = (1 - gamma_0) plus/y
		// with 1 - gamma_0:
		// | 1  0  1  0 |        | psi.e0 + psi.e2 |
		// | 0  1  0  1 |  psi = | psi.e1 + psi.e3 |
		// | 1  0  1  0 |        | psi.e1 + psi.e3 |
		// | 0  1  0  1 |        | psi.e0 + psi.e2 |
		///////////////////////////////////
		//here, only the independent components are needed
		//the other ones are incorporated in the dirac trace
		// psia = 0. component of (1-gamma_0)plus
		psia = su3vec_acc(plus.e0, plus.e2);
		// psib = 1. component of (1-gamma_0)plus
		psib = su3vec_acc(plus.e1, plus.e3);
		phia = su3vec_acc(y.e0, y.e2);
		phib = su3vec_acc(y.e1, y.e3);
		// v1 = Tr(phi*psi_dagger)
		v1 = tr_v_times_u_dagger(phia, psia, phib, psib);
		//U*v1 = U*(phi_a)
		tmp = matrix_su3to3x3(U);
		v2 = multiply_matrix3x3_dagger (tmp, v1);
		v1 = multiply_matrix3x3_by_complex(v2, bc_tmp);
		out_tmp = tr_lambda_u(v1);
		update_gaugemomentum(out_tmp, 1., global_link_pos, out);

		/////////////////////////////////////
		//mu = -0
		nn = get_lower_neighbor_temporal(t);
		global_link_pos_down = get_link_pos(dir, n, nn);
		//transform normal indices to eoprec index
		nn_eo = get_n_eoprec(n, nn);
		plus = getSpinor_eo(X, nn_eo);
		U = get_matrixsu3(field, n, nn, dir);
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
		///////////////////////////////////
		// Calculate psi/phi = (1 + gamma_0) y
		// with 1 + gamma_0:
		// | 1  0 -1  0 |       | psi.e0 - psi.e2 |
		// | 0  1  0 -1 | psi = | psi.e1 - psi.e3 |
		// |-1  0  1  0 |       | psi.e1 - psi.e2 |
		// | 0 -1  0  1 |       | psi.e0 - psi.e3 |
		///////////////////////////////////
		psia = su3vec_dim(plus.e0, plus.e2);
		psib = su3vec_dim(plus.e1, plus.e3);
		phia = su3vec_dim(y.e0, y.e2);
		phib = su3vec_dim(y.e1, y.e3);
		//CP: here is the difference with regard to +mu-direction: psi and phi interchanged!!
		// v1 = Tr(psi*phi_dagger)
		v1 = tr_v_times_u_dagger(psia, phia, psib, phib);
		//U*v1 = U*(phi_a)
		tmp = matrix_su3to3x3(U);
		v2 = multiply_matrix3x3_dagger (tmp, v1);
		v1 = multiply_matrix3x3_by_complex(v2, bc_tmp);
		out_tmp = tr_lambda_u(v1);
		update_gaugemomentum(out_tmp, 1., global_link_pos_down, out);

		//comments correspond to the mu=0 ones
		/////////////////////////////////
		//mu = 1
		/////////////////////////////////
		dir = 1;
		//this stays the same for all spatial directions at the moment
		bc_tmp.re = 2.* kappa_in * SPATIAL_RE;
		bc_tmp.im = 2.* kappa_in * SPATIAL_IM;

		/////////////////////////////////
		//mu = +1
		global_link_pos = get_link_pos(dir, n, t);
		nn = get_neighbor(n, dir);
		//transform normal indices to eoprec index
		nn_eo = get_n_eoprec(nn, t);
		plus = getSpinor_eo(X, nn_eo);
		U = get_matrixsu3(field, n, t, dir);
		/////////////////////////////////
		//Calculate (1 - gamma_1) y
		//with 1 - gamma_1:
		//| 1  0  0  i |       |       psi.e0 + i*psi.e3  |
		//| 0  1  i  0 | psi = |       psi.e1 + i*psi.e2  |
		//| 0  i  1  0 |       |(-i)*( psi.e1 + i*psi.e2) |
		//| i  0  0  1 |       |(-i)*( psi.e0 + i*psi.e3) |
		/////////////////////////////////
		//psi = (1-gamma_mu)plus
		psia = su3vec_acc_i(plus.e0, plus.e3);
		psib = su3vec_acc_i(plus.e1, plus.e2);
		phia = su3vec_acc_i(y.e0, y.e3);
		phib = su3vec_acc_i(y.e1, y.e2);

		v1 = tr_v_times_u_dagger(phia, psia, phib, psib);
		tmp = matrix_su3to3x3(U);
		v2 = multiply_matrix3x3_dagger (tmp, v1);
		v1 = multiply_matrix3x3_by_complex(v2, bc_tmp);
		out_tmp = tr_lambda_u(v1);
		update_gaugemomentum(out_tmp, 1., global_link_pos, out);
		///////////////////////////////////
		//mu = -1
		nn = get_lower_neighbor(n, dir);
		global_link_pos_down = get_link_pos(dir, nn, t);
		//transform normal indices to eoprec index
		nn_eo = get_n_eoprec(nn, t);
		plus = getSpinor_eo(X, nn_eo);
		U = get_matrixsu3(field, nn, t, dir);
		///////////////////////////////////
		// Calculate (1 + gamma_1) y
		// with 1 + gamma_1:
		// | 1  0  0 -i |       |       psi.e0 - i*psi.e3  |
		// | 0  1 -i  0 | psi = |       psi.e1 - i*psi.e2  |
		// | 0  i  1  0 |       |(-i)*( psi.e1 - i*psi.e2) |
		// | i  0  0  1 |       |(-i)*( psi.e0 - i*psi.e3) |
		///////////////////////////////////
		psia = su3vec_dim_i(plus.e0, plus.e3);
		psib = su3vec_dim_i(plus.e1, plus.e2);
		phia = su3vec_dim_i(y.e0, y.e3);
		phib = su3vec_dim_i(y.e1, y.e2);

		v1 = tr_v_times_u_dagger(psia, phia, psib, phib);
		tmp = matrix_su3to3x3(U);
		v2 = multiply_matrix3x3_dagger (tmp, v1);
		v1 = multiply_matrix3x3_by_complex(v2, bc_tmp);
		out_tmp = tr_lambda_u(v1);
		update_gaugemomentum(out_tmp, 1., global_link_pos_down, out);

		/////////////////////////////////
		//mu = 2
		/////////////////////////////////
		dir = 2;

		///////////////////////////////////
		// mu = +2
		global_link_pos = get_link_pos(dir, n, t);
		nn = get_neighbor(n, dir);
		//transform normal indices to eoprec index
		nn_eo = get_n_eoprec(nn, t);
		plus = getSpinor_eo(X, nn_eo);
		U = get_matrixsu3(field, n, t, dir);
		///////////////////////////////////
		// Calculate (1 - gamma_2) y
		// with 1 - gamma_2:
		// | 1  0  0  1 |       |       psi.e0 + psi.e3  |
		// | 0  1 -1  0 | psi = |       psi.e1 - psi.e2  |
		// | 0 -1  1  0 |       |(-1)*( psi.e1 + psi.e2) |
		// | 1  0  0  1 |       |     ( psi.e0 + psi.e3) |
		///////////////////////////////////
		psia = su3vec_acc(plus.e0, plus.e3);
		psib = su3vec_dim(plus.e1, plus.e2);
		phia = su3vec_acc(y.e0, y.e3);
		phib = su3vec_dim(y.e1, y.e2);

		v1 = tr_v_times_u_dagger(phia, psia, phib, psib);
		tmp = matrix_su3to3x3(U);
		v2 = multiply_matrix3x3_dagger (tmp, v1);
		v1 = multiply_matrix3x3_by_complex(v2, bc_tmp);
		out_tmp = tr_lambda_u(v1);
		update_gaugemomentum(out_tmp, 1., global_link_pos, out);

		///////////////////////////////////
		//mu = -2
		nn = get_lower_neighbor(n, dir);
		global_link_pos_down = get_link_pos(dir, nn, t);
		//transform normal indices to eoprec index
		nn_eo = get_n_eoprec(nn, t);
		plus = getSpinor_eo(X, nn_eo);
		U = get_matrixsu3(field, nn, t, dir);

		///////////////////////////////////
		// Calculate (1 + gamma_2) y
		// with 1 + gamma_2:
		// | 1  0  0 -1 |       |       psi.e0 - psi.e3  |
		// | 0  1  1  0 | psi = |       psi.e1 + psi.e2  |
		// | 0  1  1  0 |       |     ( psi.e1 + psi.e2) |
		// |-1  0  0  1 |       |(-1)*( psi.e0 - psi.e3) |
		///////////////////////////////////
		psia = su3vec_dim(plus.e0, plus.e3);
		psib = su3vec_acc(plus.e1, plus.e2);
		phia = su3vec_dim(y.e0, y.e3);
		phib = su3vec_acc(y.e1, y.e2);

		v1 = tr_v_times_u_dagger(psia, phia, psib, phib);
		tmp = matrix_su3to3x3(U);
		v2 = multiply_matrix3x3_dagger (tmp, v1);
		v1 = multiply_matrix3x3_by_complex(v2, bc_tmp);
		out_tmp = tr_lambda_u(v1);
		update_gaugemomentum(out_tmp, 1., global_link_pos_down, out);

		/////////////////////////////////
		//mu = 3
		/////////////////////////////////
		dir = 3;

		///////////////////////////////////
		// mu = +3
		global_link_pos = get_link_pos(dir, n, t);
		nn = get_neighbor(n, dir);
		//transform normal indices to eoprec index
		nn_eo = get_n_eoprec(nn, t);
		plus = getSpinor_eo(X, nn_eo);
		U = get_matrixsu3(field, n, t, dir);

		///////////////////////////////////
		// Calculate (1 - gamma_3) y
		// with 1 - gamma_3:
		// | 1  0  i  0 |        |       psi.e0 + i*psi.e2  |
		// | 0  1  0 -i |  psi = |       psi.e1 - i*psi.e3  |
		// |-i  0  1  0 |        |   i *(psi.e0 + i*psi.e2) |
		// | 0  i  0  1 |        | (-i)*(psi.e1 - i*psi.e3) |
		///////////////////////////////////
		psia = su3vec_acc_i(plus.e0, plus.e2);
		psib = su3vec_dim_i(plus.e1, plus.e3);
		phia = su3vec_acc_i(y.e0, y.e2);
		phib = su3vec_dim_i(y.e1, y.e3);

		v1 = tr_v_times_u_dagger(phia, psia, phib, psib);
		tmp = matrix_su3to3x3(U);
		v2 = multiply_matrix3x3_dagger (tmp, v1);
		v1 = multiply_matrix3x3_by_complex(v2, bc_tmp);
		out_tmp = tr_lambda_u(v1);
		update_gaugemomentum(out_tmp, 1., global_link_pos, out);

		///////////////////////////////////
		//mu = -3
		nn = get_lower_neighbor(n, dir);
		global_link_pos_down = get_link_pos(dir, nn, t);
		//transform normal indices to eoprec index
		nn_eo = get_n_eoprec(nn, t);
		plus = getSpinor_eo(X, nn_eo);
		U = get_matrixsu3(field, nn, t, dir);

		///////////////////////////////////
		// Calculate (1 + gamma_3) y
		// with 1 + gamma_3:
		// | 1  0 -i  0 |       |       psi.e0 - i*psi.e2  |
		// | 0  1  0  i | psi = |       psi.e1 + i*psi.e3  |
		// | i  0  1  0 |       | (-i)*(psi.e0 - i*psi.e2) |
		// | 0 -i  0  1 |       |   i *(psi.e1 + i*psi.e3) |
		///////////////////////////////////
		psia = su3vec_dim_i(plus.e0, plus.e2);
		psib = su3vec_acc_i(plus.e1, plus.e3);
		phia = su3vec_dim_i(y.e0, y.e2);
		phib = su3vec_acc_i(y.e1, y.e3);

		v1 = tr_v_times_u_dagger(psia, phia, psib, phib);
		tmp = matrix_su3to3x3(U);
		v2 = multiply_matrix3x3_dagger (tmp, v1);
		v1 = multiply_matrix3x3_by_complex(v2, bc_tmp);
		out_tmp = tr_lambda_u(v1);
		update_gaugemomentum(out_tmp, 1., global_link_pos_down, out);
	}
}
