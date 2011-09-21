__kernel void M_tm_sitediagonal(__global spinorfield_eoprec * in, __global spinorfield_eoprec * out){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	spinor out_tmp;
	spinor plus;
	
	hmc_complex twistfactor = {1., MUBAR};
	hmc_complex twistfactor_minus = {1., MMUBAR};
	
	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE; id_tmp += global_size) {	
		out_tmp = set_spinor_zero();
		//get input spinor
		plus = get_spinor_from_eoprec_field(in, id_tmp);
		//Diagonalpart:
		out_tmp = M_diag_tm_local(plus, twistfactor, twistfactor_minus);
		put_spinor_to_eoprec_field(out_tmp, out, id_tmp);
	}
}

__kernel void M_tm_inverse_sitediagonal(__global spinorfield_eoprec * in, __global spinorfield_eoprec * out){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	spinor out_tmp;
	spinor plus;
	hmc_complex twistfactor = {1., MUBAR};
	hmc_complex twistfactor_minus = {1., MMUBAR};
	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE; id_tmp += global_size) {	
		out_tmp = set_spinor_zero();
		//get input spinor
		plus = get_spinor_from_eoprec_field(in, id_tmp);
		//Diagonalpart, here the twisted factor give the inverse matrix:
		out_tmp = M_diag_tm_local(plus, twistfactor_minus, twistfactor);
		hmc_float denom = 1./(1. + MUBAR*MUBAR);
		out_tmp = real_multiply_spinor(out_tmp,denom);
		put_spinor_to_eoprec_field(out_tmp, out, id_tmp);
	}
}

__kernel void M_tm_sitediagonal_minus(__global spinorfield_eoprec * in, __global spinorfield_eoprec * out){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	spinor out_tmp;
	spinor plus;
	
	hmc_complex twistfactor = {1., MUBAR};
	hmc_complex twistfactor_minus = {1., MMUBAR};
	
	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE; id_tmp += global_size) {	
		out_tmp = set_spinor_zero();
		//get input spinor
		plus = get_spinor_from_eoprec_field(in, id_tmp);
		//Diagonalpart:
		out_tmp = M_diag_tm_local(plus, twistfactor_minus, twistfactor);
		put_spinor_to_eoprec_field(out_tmp, out, id_tmp);
	}
}

__kernel void M_tm_inverse_sitediagonal_minus(__global spinorfield_eoprec * in, __global spinorfield_eoprec * out){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	spinor out_tmp;
	spinor plus;
	hmc_complex twistfactor = {1., MUBAR};
	hmc_complex twistfactor_minus = {1., MMUBAR};
	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE; id_tmp += global_size) {	
		out_tmp = set_spinor_zero();
		//get input spinor
		plus = get_spinor_from_eoprec_field(in, id_tmp);
		//Diagonalpart, here the twisted factor give the inverse matrix:
		out_tmp = M_diag_tm_local(plus, twistfactor, twistfactor_minus);
		hmc_float denom = 1./(1. + MUBAR*MUBAR);
		out_tmp = real_multiply_spinor(out_tmp,denom);
		put_spinor_to_eoprec_field(out_tmp, out, id_tmp);
	}
}

