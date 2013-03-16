/**
 * @file operations used for the molecular dynamics update
 */
/**
 * @file operations used for the molecular dynamics update of the gaugefield
 * u_out = exp(i eps p_in) u_in
 */

__kernel void md_update_gaugefield(const hmc_float eps, __global const aeStorageType * const restrict p_in, __global Matrixsu3StorageType * const restrict u_inout)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int index;
	Matrixsu3 tmp;
	Matrixsu3 tmp2;

	//CP: it is GAUGEMOMENTASIZE = NDIM * SPINORFIELDSIZE
	for(int id_local = id; id_local < VOL4D_LOCAL; id_local += global_size) {
		/** @todo this must be done more efficient */
		st_index pos = (id_local < VOL4D_LOCAL / 2) ? get_even_st_idx_local(id_local) : get_odd_st_idx_local(id_local - (VOL4D_LOCAL / 2));
		for(int mu = 0; mu < NDIM; mu++) {
			index = get_link_idx(mu, pos);
			// an su3 algebra element has NC*NC-1 = 8 hmc_float entries
			// &(p_in[index*8]) should point to the right position for the pos-th element of the long gaugemomentum vector p_in
			tmp2 = build_su3matrix_by_exponentiation(getAe(p_in, index), eps);
			tmp2 = project_su3(tmp2);
			tmp = get_matrixsu3(u_inout, pos.space, pos.time, mu);
			tmp2 = multiply_matrixsu3( tmp2, tmp);
			put_matrixsu3(u_inout, tmp2, pos.space, pos.time, mu);
		}
	}
}
