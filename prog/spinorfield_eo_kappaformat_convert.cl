__kernel void convert_to_kappa_format_eoprec( __global spinorfield_eoprec *in)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	hmc_float tmp = sqrt(KAPPA*2.);

	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE2; id_tmp += global_size) {
		in[id_tmp] = real_multiply_spinor(in[id_tmp], tmp);
	}
	return;
}

__kernel void convert_from_kappa_format_eoprec( __global spinorfield_eoprec *in, __global spinorfield_eoprec* out)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	hmc_float tmp = 1./sqrt(KAPPA*2.);

	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE2; id_tmp += global_size) {
		in[id_tmp] = real_multiply_spinor(in[id_tmp], tmp);
	}
	return;
}
