__kernel void localQ_test(__global ocl_s_gaugefield * field, __global hmc_float * out)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//CP: this is essentially the plaquette-kernel. The result is however different since the normalization is missing
	for(int id_tmp = id; id_tmp < VOL4D; id_tmp += global_size) {
		st_index pos = (id_tmp % 2 == 0) ? get_even_site(id_tmp / 2) : get_odd_site(id_tmp / 2);
		Matrix3x3 V;
		hmc_float res = 0.;

		//NOTE: The kernel crashes on ATI GPUs if one start with mu=0 here, although this does not make any sense since the nu-loop then does not contribute!!!
		for(int mu = 1; mu < NDIM; mu++) {
			for(int nu = 0; nu < mu; nu++) {
				Matrix3x3 tmp;
				tmp = 	local_Q_plaquette(field,pos.space,pos.time, mu, nu);
				res += trace_matrix3x3(tmp).re;
		}}

		int global_pos = get_global_pos(pos.space, pos.time);
		out[global_pos] = res;
	}

}
