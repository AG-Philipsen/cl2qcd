/**
 * @file operations used for the molecular dynamics update
 */

//CP: molecular dynamics update for the gauge momenta:
//p_out = p_in - eps/2 force(u_in, phi)
//it is assumed that the force term has already been computed. then one only has real-vectors and this is essentially adding one vector to another...
__kernel void md_update_gaugemomenta(hmc_float eps, __global ae * p_inout, __global ae* force_in)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	for(int id_tmp = id; id_tmp < GAUGEMOMENTASIZE; id_tmp += global_size) {
		update_gaugemomentum(force_in[id_tmp], eps, id_tmp, p_inout);
	}

}
