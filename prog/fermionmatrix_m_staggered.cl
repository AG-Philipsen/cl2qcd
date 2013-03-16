/**
 * @file M normal staggered fermionmatrix
 */


__kernel void M_staggered(__global const su3vec * const restrict in, __global const Matrixsu3StorageType * const restrict field, __global su3vec * const restrict out, hmc_float mass_in)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	int n, t;
	su3vec out_tmp, out_tmp2;

	for(int id_local = id; id_local < SPINORFIELDSIZE_LOCAL; id_local += global_size) {

		/** @todo this must be done more efficient */
		st_index pos = (id_local % 2 == 0) ? get_even_st_idx_local(id_local / 2) : get_odd_st_idx_local(id_local / 2);

		//this assumes the form: M = m - dslash

		//Diagonalpart: m * chi(n)
		out_tmp = get_su3vec_from_field(in, pos.space, pos.time);
		out_tmp = su3vec_times_real(out_tmp, mass_in);
		//calc dslash
		out_tmp2 = dslash_local_0(in, field, pos.space, pos.time);
		out_tmp = su3vec_dim(out_tmp, out_tmp2);
		// these are not yet implemented...
		/*
		out_tmp2 = dslash_local_1(in, field, pos.space, pos.time);
		out_tmp = su3vec_dim(out_tmp, out_tmp2);
		out_tmp2 = dslash_local_2(in, field, pos.space, pos.time);
		out_tmp = su3vec_dim(out_tmp, out_tmp2);
		out_tmp2 = dslash_local_3(in, field, pos.space, pos.time);
		out_tmp = su3vec_dim(out_tmp, out_tmp2);
		*/
		
		put_su3vec_to_field(out_tmp, out, pos.space, pos.time);
	}
}
