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
	Matrixsu3 tmp;
	Matrixsu3 tmp2;

	//CP: it is GAUGEMOMENTASIZE = NDIM * SPINORFIELDSIZE
	PARALLEL_FOR(index, VOL4D_MEM * NDIM) {
		// an su3 algebra element has NC*NC-1 = 8 hmc_float entries
		// &(p_in[index*8]) should point to the right position for the pos-th element of the long gaugemomentum vector p_in
		tmp2 = build_su3matrix_by_exponentiation(getAe(p_in, index), eps);
		tmp2 = project_su3(tmp2);
		tmp = getSU3(u_inout, index);
		tmp2 = multiply_matrixsu3(tmp2, tmp);
		putSU3(u_inout, index, tmp2);
	}
}
