__kernel void set_zero_spinorfield( __global spinor * const restrict x )
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size) {
		x[id_tmp] = set_spinor_zero();
	}
	return;
}
