// Description of variables of saxpy:
//  - x: The first input staggered field (an su3vec per each site => vector of VOL4D/2
//       components that are su3vec varibles)
//  - y: The second input staggered field (an su3vec per each site => vector of VOL4D/2
//       components that are su3vec varibles)
//  - alpha: The complex number by which x has to be multiplied
//  - out: The output staggered field: alpha*x+y (site by site)

__kernel void saxpy_staggered_eoprec(__global const staggeredStorageType * const x, __global const staggeredStorageType * const y, __global const hmc_complex * const alpha, __global staggeredStorageType * const out)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	const hmc_complex alpha_tmp = complexLoadHack(alpha);
	for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {
		su3vec x_tmp = get_su3vec_from_field_eo(x, id_mem);
		su3vec y_tmp = get_su3vec_from_field_eo(y, id_mem);
		x_tmp = su3vec_times_complex(x_tmp, alpha_tmp);
		x_tmp = su3vec_acc(y_tmp, x_tmp);
		put_su3vec_to_field_eo(out, id_mem, x_tmp);
	}
}

__kernel void saxpy_arg_staggered_eoprec(__global const spinorStorageType * const x, __global const spinorStorageType * const y, const hmc_float alpha_re, const hmc_float alpha_im, __global spinorStorageType * const out)
{
	const int id = get_global_id(0);
	const int global_size = get_global_size(0);

	const hmc_complex alpha = (hmc_complex) {alpha_re, alpha_im};

	for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {
		su3vec x_tmp = get_su3vec_from_field_eo(x, id_mem);
		su3vec y_tmp = get_su3vec_from_field_eo(y, id_mem);
		x_tmp = su3vec_times_complex(x_tmp, alpha);
		x_tmp = su3vec_acc(y_tmp, x_tmp);
		put_su3vec_to_field_eo(out, id_mem, x_tmp);
	}
}
