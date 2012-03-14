/**
 @file fermionmatrix-functions for eoprec spinorfields
*/

//"local" dslash working on a particular link (n,t) of an eoprec field
//NOTE: each component is multiplied by +KAPPA, so the resulting spinor has to be mutliplied by -1 to obtain the correct dslash!!!
//the difference to the "normal" dslash is that the coordinates of the neighbors have to be transformed into an eoprec index
spinor dslash_eoprec_unified_local(__global const spinorStorageType * const restrict in, __global Matrixsu3StorageType  const * const restrict field, const st_idx idx_arg, const dir_idx dir, hmc_float kappa_in)
{
	//this is used to save the idx of the neighbors
	st_idx idx_neigh;

	spinor out_tmp, plus;
	site_idx nn_eo;
	su3vec psi, phi;
	Matrixsu3 U;
	//this is used to save the BC-conditions...
	hmc_complex bc_tmp = (dir == TDIR) ? (hmc_complex) {
		kappa_in * TEMPORAL_RE, kappa_in * TEMPORAL_IM
} :
	(hmc_complex) {
		kappa_in * SPATIAL_RE, kappa_in * SPATIAL_IM
	};
	out_tmp = set_spinor_zero();

	///////////////////////////////////
	// mu = +dir
	idx_neigh = get_neighbor_from_st_idx(idx_arg, dir);
	//transform normal indices to eoprec index
	nn_eo = get_eo_site_idx_from_st_idx(idx_neigh);
	plus = getSpinor_eo(in, nn_eo);
	U = getSU3(field, get_link_idx_SOA(dir, idx_arg));
	if(dir == XDIR) {
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
	} else if(dir == YDIR) {
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
	} else if(dir == ZDIR) {
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
	} else { // TDIR
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
	}

	///////////////////////////////////
	//mu = -dir
	idx_neigh = get_lower_neighbor_from_st_idx(idx_arg, dir);
	//transform normal indices to eoprec index
	nn_eo = get_eo_site_idx_from_st_idx(idx_neigh);
	plus = getSpinor_eo(in, nn_eo);
	U = getSU3(field, get_link_idx_SOA(dir, idx_neigh));
	//in direction -mu, one has to take the complex-conjugated value of bc_tmp. this is done right here.
	bc_tmp = (dir == TDIR) ? (hmc_complex) {
		kappa_in * TEMPORAL_RE, 1. * kappa_in * TEMPORAL_IM
} :
	(hmc_complex) {
		kappa_in * SPATIAL_RE, -1.*kappa_in * SPATIAL_IM
	};
	if(dir == XDIR) {
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
	} else if(dir == YDIR) {
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
	} else if(dir == ZDIR) {
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
	} else { // TDIR
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
		///////////////////////////////////
		// Calculate psi/phi = (1 + gamma_0) y
		// with 1 + gamma_0:
		// | 1  0 -1  0 |       | psi.e0 - psi.e2 |
		// | 0  1  0 -1 | psi = | psi.e1 - psi.e3 |
		// |-1  0  1  0 |       | psi.e1 - psi.e2 |
		// | 0 -1  0  1 |       | psi.e0 - psi.e3 |
		///////////////////////////////////
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
	}

	return out_tmp;
}

