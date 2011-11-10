/**
 * @file stout-smearing-fermion-force
 */

//this is done after the tmlqcd-analogue (stout_smear_force.c)
__kernel void stout_smear_fermion_force()
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);


}
