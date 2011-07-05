__kernel void gamma5_eoprec(__global spinorfield_eoprec *in, __global spinorfield_eoprec *out){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	spinor out_tmp;

	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE; id_tmp += global_size) {	

		out_tmp = get_spinor_from_eoprec_field(in, id_tmp); 
		out_tmp = gamma5_local(out_tmp);
		put_spinor_to_eoprec_field(out_tmp, out, id_tmp);
	}
}
