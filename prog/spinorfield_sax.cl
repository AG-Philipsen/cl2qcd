// out = alpha*x
__kernel void sax(__global const spinor * const restrict x, __global const hmc_complex * const restrict alpha, __global spinor * const restrict out)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	hmc_complex alpha_tmp = complexLoadHack(alpha);
	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size) {
		spinor x_tmp = x[id_tmp];
		out[id_tmp] = spinor_times_complex(x_tmp, alpha_tmp);
	}
}
