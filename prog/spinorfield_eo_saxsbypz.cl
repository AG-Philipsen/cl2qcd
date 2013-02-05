//alpha*x + beta*y + z
__kernel void saxsbypz_eoprec(__global const spinorStorageType * const x, __global const spinorStorageType * const y, __global const spinorStorageType * const z, __global const hmc_complex * const alpha, __global hmc_complex * beta, __global spinorStorageType * const out)
{
	const int id = get_global_id(0);
	const int global_size = get_global_size(0);

	const hmc_complex alpha_tmp = complexLoadHack(alpha);
	const hmc_complex beta_tmp = complexLoadHack(beta);
	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE; id_tmp += global_size) {
		spinor x_tmp = getSpinor_eo(x, id_tmp);
		x_tmp = spinor_times_complex(x_tmp, alpha_tmp);
		spinor y_tmp = getSpinor_eo(y, id_tmp);
		y_tmp = spinor_times_complex(y_tmp, beta_tmp);
		spinor z_tmp = getSpinor_eo(z, id_tmp);

		spinor out_tmp = spinor_acc_acc(y_tmp, x_tmp, z_tmp);
		putSpinor_eo(out, id_tmp, out_tmp);
	}

	return;
}
