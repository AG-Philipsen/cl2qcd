//this calculates the sum of all elements in the matrix
hmc_float sum_up_matrix3x3(Matrix3x3 in)
{
	return in.e00.re + in.e00.im +
	       in.e01.re + in.e01.im +
	       in.e02.re + in.e02.im +
	       in.e10.re + in.e10.im +
	       in.e11.re + in.e11.im +
	       in.e12.re + in.e12.im +
	       in.e20.re + in.e20.im +
	       in.e21.re + in.e21.im +
	       in.e22.re + in.e22.im ;
}


__kernel void staple_test(__global Matrixsu3StorageType * field, __global hmc_float * out)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);

	for(int id_tmp = id; id_tmp < VOL4D; id_tmp += global_size) {
		st_index pos = (id_tmp % 2 == 0) ? get_even_site(id_tmp / 2) : get_odd_site(id_tmp / 2);
		Matrix3x3 V;
		hmc_float res = 0;
		hmc_float res2;
		for(int mu = 0; mu < NDIM; mu++) {
			V = calc_staple(field, pos.space, pos.time, mu);
			res2 = sum_up_matrix3x3(V);
			res += res2;
		}
		int global_pos = get_global_pos(pos.space, pos.time);
		out[global_pos] = res;
	}

}
