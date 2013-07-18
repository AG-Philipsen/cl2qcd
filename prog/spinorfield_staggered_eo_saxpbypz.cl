// Description of variables of saxpbypz:
//  - x: The first input staggered field (an su3vec per each site => vector of VOL4D/2
//       components that are su3vec varibles)
//  - y: The second input staggered field (an su3vec per each site => vector of VOL4D/2
//       components that are su3vec varibles)
//  - z: The third input staggered field (an su3vec per each site => vector of VOL4D/2
//       components that are su3vec varibles)
//  - alpha: The complex number by which x has to be multiplied
//  - beta:  The complex number by which y has to be multiplied
//  - out: The output staggered field: alpha*x+beta*y+z (site by site)

__kernel void saxpbypz_staggered_eoprec(__global const staggeredStorageType * const x, __global const staggeredStorageType * const y, __global const staggeredStorageType * const z, __global const hmc_complex * const alpha, __global hmc_complex * beta, __global staggeredStorageType * const out)
{
	const int id = get_global_id(0);
	const int global_size = get_global_size(0);

	const hmc_complex alpha_tmp = complexLoadHack(alpha);
	const hmc_complex beta_tmp = complexLoadHack(beta);
	for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {
		su3vec x_tmp = get_su3vec_from_field_eo(x, id_mem);
		x_tmp = su3vec_times_complex(x_tmp, alpha_tmp);
		su3vec y_tmp = get_su3vec_from_field_eo(y, id_mem);
		y_tmp = su3vec_times_complex(y_tmp, beta_tmp);
		su3vec z_tmp = get_su3vec_from_field_eo(z, id_mem);

		su3vec out_tmp = su3vec_acc_acc(y_tmp, x_tmp, z_tmp);
		put_su3vec_to_field_eo(out, id_mem, out_tmp);
	}

	return;
}
