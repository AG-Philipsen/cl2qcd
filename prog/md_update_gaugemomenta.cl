/**
 * @file operations used for the molecular dynamics update of the gauge momenta
 * p_out = p_in - eps/2 force(u_in, phi)
 * It is assumed that the force term has already been computed. Then one only has real-vectors and this is essentially adding one vector to another...
 */

__kernel void md_update_gaugemomenta(const hmc_float eps, __global aeStorageType * const restrict p_inout, __global const aeStorageType * const restrict force_in)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	for(int id_tmp = id; id_tmp < GAUGEMOMENTASIZE; id_tmp += global_size) {
		update_gaugemomentum(getAe(force_in, id_tmp), eps, id_tmp, p_inout);
	}
}
