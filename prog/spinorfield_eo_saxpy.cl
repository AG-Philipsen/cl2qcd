// -alpha*x + y
//CP: defined with a minus!!!
__kernel void saxpy_eoprec(__global const spinorStorageType * const x, __global const spinorStorageType * const y, __global const hmc_complex * const alpha, __global spinorStorageType * const out)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	const hmc_complex alpha_tmp = complexLoadHack(alpha);
	for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {
		spinor x_tmp = getSpinor_eo(x, id_mem);
		spinor y_tmp = getSpinor_eo(y, id_mem);
		x_tmp = spinor_times_complex(x_tmp, alpha_tmp);
		x_tmp = spinor_dim(y_tmp, x_tmp);
		putSpinor_eo(out, id_mem, x_tmp);
	}
}

__kernel void saxpy_arg_eoprec(__global const spinorStorageType * const x, __global const spinorStorageType * const y, const hmc_float alpha_re, const hmc_float alpha_im, __global spinorStorageType * const out)
{
	const int id = get_global_id(0);
	const int global_size = get_global_size(0);

	const hmc_complex alpha = (hmc_complex) {
		alpha_re, alpha_im
	};

	for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {
		spinor x_tmp = getSpinor_eo(x, id_mem);
		spinor y_tmp = getSpinor_eo(y, id_mem);
		x_tmp = spinor_times_complex(x_tmp, alpha);
		x_tmp = spinor_dim(y_tmp, x_tmp);
		putSpinor_eo(out, id_mem, x_tmp);
	}
}
