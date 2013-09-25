/**
 * @file
 * Kernel for the gaugemomentum saxpy operation.
 *  @param x The first input gaugemomentum field (an ae per each site)
 *  @param y The second input gaugemomentum field (an su3vec per each site)
 *  @param alpha The REAL number which x has to be multiplied by
 *  @param out The output gaugemomentum field: alpha*x+y (site by site)
 */
__kernel void gaugemomentum_saxpy(__global const aeStorageType * const x, __global const aeStorageType * const y, __global const hmc_float * const alpha, __global aeStorageType * const out)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	for(int id_mem = id; id_mem < GAUGEMOMENTASIZE_MEM; id_mem += global_size) {
		ae x_tmp = getAe(x, id_mem);
		ae y_tmp = getAe(y, id_mem);
		putAe(out, id_mem, acc_factor_times_algebraelement(y_tmp, *alpha, x_tmp));
	}
}



#if 0
/**
 * @file operations used for the molecular dynamics update of the gauge momenta
 * p_out = p_in - eps/2 force(u_in, phi)
 * It is assumed that the force term has already been computed. Then one only has real-vectors and this is essentially adding one vector to another...
 */
__kernel void md_update_gaugemomenta(const hmc_float eps, __global aeStorageType * const restrict p_inout, __global const aeStorageType * const restrict force_in)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	for(int id_mem = id; id_mem < GAUGEMOMENTASIZE_MEM; id_mem += global_size) {
		update_gaugemomentum(getAe(force_in, id_mem), eps, id_mem, p_inout);
	}
}
#endif