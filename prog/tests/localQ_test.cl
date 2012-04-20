__kernel void localQ_test(__global Matrixsu3StorageType * field, __global hmc_float * out)
{
	//CP: this is essentially the plaquette-kernel. The result is however different since the normalization is missing
	PARALLEL_FOR(id_tmp, VOL4D) {
		st_index pos = (id_tmp % 2 == 0) ? get_even_site(id_tmp / 2) : get_odd_site(id_tmp / 2);
		hmc_float res = 0.;

		//NOTE: The kernel crashes on ATI GPUs if one start with mu=0 here, although this does not make any sense since the nu-loop then does not contribute!!!
		for(int mu = 1; mu < NDIM; mu++) {
			for(int nu = 0; nu < mu; nu++) {
				Matrix3x3 tmp;
				tmp = local_Q_plaquette(field, pos.space, pos.time, mu, nu);
				res += trace_matrix3x3(tmp).re;
			}
		}

		int global_pos = get_global_pos(pos.space, pos.time);
		out[global_pos] = res;
	}

}
