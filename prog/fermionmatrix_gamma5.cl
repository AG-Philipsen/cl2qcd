__kernel void gamma5(__global spinorfield *in, __global spinorfield *out){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	spinor out_tmp;
	int n,t;

	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size) {	
		/** @todo this must be done more efficient */
		if(id_tmp%2 == 0) get_even_site(id_tmp/2, &n, &t);
		else get_odd_site(id_tmp/2, &n, &t);

		out_tmp = get_spinor_from_field(in,n, t); 
		out_tmp = gamma5_local(out_tmp);
		put_spinor_to_field(out_tmp, out, n, t);
	}
}
