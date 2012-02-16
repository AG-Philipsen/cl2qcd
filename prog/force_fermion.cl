/**
 * @file kernel for the non-eo fermion force
 */

__kernel void fermion_force(__global ocl_s_gaugefield * field, __global  spinorfield * Y, __global  spinorfield * X, __global  ae * out)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	for(int global_link_pos = id; global_link_pos < GAUGEMOMENTASIZE; global_link_pos += global_size) {
		int dir = global_link_pos / SPINORFIELDSIZE;
		int site_idx = global_link_pos % SPINORFIELDSIZE;

		st_index pos = (site_idx % 2 == 0) ? get_even_site(site_idx / 2) : get_odd_site(site_idx / 2);

		Matrixsu3 U;
		Matrix3x3 v1, v2, tmp;
		su3vec psia, psib, phia, phib;
		spinor y, plus;
		int nn;
		ae out_tmp;
		//this is used to save the BC-conditions...
		hmc_complex bc_tmp;
		int n = pos.space;
		int t = pos.time;

		if(dir == 0) {
			y = get_spinor_from_field(Y, n, t);
			///////////////////////////////////
			// Calculate gamma_5 y
			///////////////////////////////////
			y = gamma5_local(y);

			///////////////////////////////////
			// mu = 0
			///////////////////////////////////
			nn = get_neighbor_temporal(t);

			//the 2 here comes from Tr(lambda_ij) = 2delta_ij
			bc_tmp.re = 2.*KAPPA_TEMPORAL_RE;
			bc_tmp.im = 2.*KAPPA_TEMPORAL_IM;

			///////////////////////////////////
			//mu = +0
			plus = get_spinor_from_field(X, n, nn);
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
			y = get_spinor_from_field(Y, n, nn);
			y = gamma5_local(y);
			plus = get_spinor_from_field(X, n, t);
			U = get_matrixsu3(field, n, t, dir);
			//if chemical potential is activated, U has to be multiplied by appropiate factor
			//this is the same as at mu=0 in the imag. case, since U is taken to be U^+ later:
			//  (exp(iq)U)^+ = exp(-iq)U^+
			//as it should be
			//in the real case, one has to take exp(q) -> exp(-q)
#ifdef _  CP_REAL_
			U = multiply_matrixsu3_by_real (U, MEXPCPR);
#endif
#ifdef _  CP_IMAG_
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
			update_gaugemomentum(out_tmp, 1., global_link_pos, out);

		} else {
			y = get_spinor_from_field(Y, n, t);
			///////////////////////////////////
			// Calculate gamma_5 y
			///////////////////////////////////
			y = gamma5_local(y);

			//comments correspond to the mu=0 ones
			/////////////////////////////////
			//mu = 1
			/////////////////////////////////
			//this stays the same for all spatial directions at the moment
			bc_tmp.re = 2.*KAPPA_SPATIAL_RE;
			bc_tmp.im = 2.*KAPPA_SPATIAL_IM;

			/////////////////////////////////
			//mu = +1
			nn = get_neighbor(n, dir);
			plus = get_spinor_from_field(X, nn, t);
			U = get_matrixsu3(field, n, t, dir);

			if( dir == 1 ) {
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
			} else if( dir == 2 ) {
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
			} else { // dir == 3
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
			}

			v1 = tr_v_times_u_dagger(phia, psia, phib, psib);
			tmp = matrix_su3to3x3(U);
			v2 = multiply_matrix3x3_dagger (tmp, v1);
			v1 = multiply_matrix3x3_by_complex(v2, bc_tmp);
			out_tmp = tr_lambda_u(v1);

			update_gaugemomentum(out_tmp, 1., global_link_pos, out);

			///////////////////////////////////
			//mu = -1
			y = get_spinor_from_field(Y, nn, t);
			y = gamma5_local(y);
			plus = get_spinor_from_field(X, n, t);
			U = get_matrixsu3(field, n, t, dir);

			if( dir == 1 ) {
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
			} else if( dir == 2 ) {
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
			} else { // dir == 3
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
			}

			v1 = tr_v_times_u_dagger(psia, phia, psib, phib);
			tmp = matrix_su3to3x3(U);
			v2 = multiply_matrix3x3_dagger (tmp, v1);
			v1 = multiply_matrix3x3_by_complex(v2, bc_tmp);

			out_tmp = tr_lambda_u(v1);
			update_gaugemomentum(out_tmp, 1., global_link_pos, out);
		}
	}
}
