/**
 * @file M normal Wilson fermionmatrix
 */

__kernel void M(__global spinorfield * in, __global ocl_s_gaugefield * field, __global spinorfield * out){
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	int n,t;
	spinor out_tmp, out_tmp2;

	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size) {	
	
		/** @todo this must be done more efficient */
		if(id_tmp%2 == 0) get_even_site(id_tmp/2, &n, &t);
		else get_odd_site(id_tmp/2, &n, &t);
		
		//Diagonalpart: (this is simple here)
		out_tmp = get_spinor_from_field(in, n, t);
		//calc dslash (this includes mutliplication with kappa)
		out_tmp2 = dslash_local(in, field, n, t);
		//M = M_diag - dslash
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		put_spinor_to_field(out_tmp, out, n, t);
	}
}
