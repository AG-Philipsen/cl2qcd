// out = alpha*x
__kernel void sax(__global spinorfield* x, __global hmc_complex * alpha, __global spinorfield* out)
{
	int id = get_global_id(0);
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id(0);

	hmc_complex alpha_tmp = (*alpha);
	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size) {
		spinor x_tmp = x[id_tmp];
		out[id_tmp] = spinor_times_complex(x_tmp, alpha_tmp);
	}
}
