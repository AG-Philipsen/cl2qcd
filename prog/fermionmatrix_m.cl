/**
 * @file M normal Wilson fermionmatrix
 */

__kernel void M_wilson(__global const spinor * const restrict in, __global const Matrixsu3StorageType * const restrict field, __global spinor * const restrict out, hmc_float kappa_in)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	spinor out_tmp, out_tmp2;

	for(int id_local = id; id_local < SPINORFIELDSIZE_LOCAL; id_local += global_size) {

		/** @todo this must be done more efficient */
		st_index pos = (id_local % 2 == 0) ? get_even_st_idx_local(id_local / 2) : get_odd_st_idx_local(id_local / 2);

		//Diagonalpart: (this is simple here)
		out_tmp = get_spinor_from_field(in, pos.space, pos.time);
		//calc dslash (this includes mutliplication with kappa)
		out_tmp2 = dslash_local_0(in, field, pos.space, pos.time, kappa_in);
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		out_tmp2 = dslash_local_1(in, field, pos.space, pos.time, kappa_in);
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		out_tmp2 = dslash_local_2(in, field, pos.space, pos.time, kappa_in);
		out_tmp = spinor_dim(out_tmp, out_tmp2);
		out_tmp2 = dslash_local_3(in, field, pos.space, pos.time, kappa_in);
		out_tmp = spinor_dim(out_tmp, out_tmp2);

		put_spinor_to_field(out_tmp, out, pos.space, pos.time);
	}
}
