__kernel void gauge_force(__global ocl_s_gaugefield * field, __global ae * out)
{

	int id = get_global_id(0);
	if(id == 0) {


		hmc_float factor = -BETA / 3.;
		int global_link_pos;
		Matrix3x3 V, tmp, tmp2;
		Matrixsu3 U;
		ae out_tmp;

		//Gauge force is factor*Im(i Tr(T_i U V))
		//   with T_i being the SU3-Generator in i-th direction and V the staplematrix
		//   and the factor being -beta/3. (for standard Wilson-action)
		for(int t = 0; t < NTIME; t++) {
			for(int n = 0; n < VOLSPACE; n++) {
				for(int mu = 0; mu < NDIM; mu++) {
					global_link_pos = get_global_link_pos(mu, n, t);
					V = calc_staple(field, n, t, mu);
					U = get_matrixsu3(field, n, t, mu);

					tmp2 = matrix_su3to3x3(U);
					tmp = multiply_matrix3x3 (tmp2, V);

					out_tmp = tr_lambda_u(tmp);

					update_gaugemomentum(out_tmp, factor, global_link_pos, out);
				}
			}
		}

	}
}
