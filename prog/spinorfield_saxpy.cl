// -alpha*x + y
//CP: defined with a minus!!!

__kernel void saxpy(__global spinor* x, __global spinor* y, __global const hmc_complex * alpha_p, __global spinor* out)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	const hmc_complex alpha = complexLoadHack(alpha_p);

	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size) {
		const spinor x_tmp = x[id_tmp];
		const spinor y_tmp = y[id_tmp];
		const spinor tmp = spinor_times_complex(x_tmp, alpha);
		out[id_tmp] = spinor_dim(y_tmp, tmp);
	}
}

// the arguments have been hacked to work on apple
__kernel void saxpy_arg(__global spinor* x, __global spinor* y, const hmc_float alpha_re, const hmc_float alpha_im, __global spinor* out)
{
	const int id = get_global_id(0);
	const int global_size = get_global_size(0);

	const hmc_complex alpha = (hmc_complex) {
		alpha_re, alpha_im
	};

	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size) {
		const spinor x_tmp = x[id_tmp];
		const spinor y_tmp = y[id_tmp];
		const spinor tmp = spinor_times_complex(x_tmp, alpha);
		out[id_tmp] = spinor_dim(y_tmp, tmp);
	}
}
