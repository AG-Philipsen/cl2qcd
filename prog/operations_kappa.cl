/** @file
 * Necessary operations regarding the measurement of the TK kappa including
 * diagonal elements of the energy-momentum-tensor
 */

//operations_kappa.cl

__kernel void T_clover_diag_partial(__global ocl_s_gaugefield* gaugefield, __global hmc_float * v_tau, __global hmc_float * v_spatial,
                                    __global hmc_float * w_tau,  __global hmc_float * w_spatial,
                                    const int mom_dir)
{

	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//Determine which diagonal elements are to be computed
	// mom_dir=3: theta_11, theta_22
	// mom_dir=2: theta_11, theta_33
	// mom_dir=1: theta_22, theta_33
	switch(mom_dir) {
		case 3:



			for(int id_tmp = id; id_tmp < VOL4D / NS; id_tmp += global_size) {

				//chose pos only on 3d subset, that means excluce z-direction
				st_index pos = (id_tmp % 2 == 0) ? get_even_site(id_tmp / 2) : get_odd_site(id_tmp / 2);

				Matrix3x3 tmp;
				tmp = local_Q_plaquette(field, pos.space, pos.time, 2, 0);
				tmp = tmp * subtract_matrix3x3_dagger(tmp, tmp);
				hmc_float tmp_20 = (trace_matrix3x3 (tmp)).re;

				tmp = local_Q_plaquette(field, pos.space, pos.time, 3, 0);
				tmp = tmp * subtract_matrix3x3_dagger(tmp, tmp);
				hmc_float tmp_30 = (trace_matrix3x3 (tmp)).re;

				tmp = local_Q_plaquette(field, pos.space, pos.time, 1, 0);
				tmp = tmp * subtract_matrix3x3_dagger(tmp, tmp);
				hmc_float tmp_10 = (trace_matrix3x3 (tmp)).re;

				//reduction
				v_tau[] += tmp_20 + tmp_30 - tmp_10;
				w_tau[] += tmp_10 + tmp_30 - tmp_20;


			}

		case 2:

		case 1:

	}

}

