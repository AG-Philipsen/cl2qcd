//alpha*x + beta*y + z
__kernel void saxsbypz(__global spinorfield* x, __global spinorfield* y, __global spinorfield* z, __global hmc_complex * alpha, __global hmc_complex * beta, __global spinorfield* out)
{
	int id = get_global_id(0);
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id(0);

	hmc_complex alpha_tmp = (*alpha);
	hmc_complex beta_tmp = (*beta);
	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size) {
		spinor x_tmp = x[id_tmp];
		x_tmp = spinor_times_complex(x_tmp, alpha_tmp);
		spinor y_tmp = y[id_tmp];
		y_tmp = spinor_times_complex(y_tmp, beta_tmp);
		spinor z_tmp = z[id_tmp];

		out[id_tmp] = spinor_acc_acc(y_tmp, x_tmp, z_tmp);
	}

	return;
}
