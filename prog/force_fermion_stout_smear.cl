/**
 * @file stout-smearing-fermion-force
 */

//this is done after the tmlqcd-analogue (stout_smear_force.c)
//I did not check this yet, but I think since the calculation of the stout-smeared force is iterative, one should most likely split this into 2 kernels...


//this kernel is meant to be the one that updates the force with all needed ingredients already computed
//	it takes the force field (out) and <fill in>
__kernel void stout_smear_fermion_force(__global  ae * out)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);


}
