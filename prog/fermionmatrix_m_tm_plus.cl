/**
 * @file M_tm_plus this is the "total" twisted-mass fermionmatrix (no evenodd) for the upper flavor
 */
__kernel void M_tm_plus(__global const spinor * const restrict in, __global const Matrixsu3StorageType * field, __global spinor * const restrict out, hmc_float kappa_in, hmc_float mubar_in)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	spinor out_tmp, out_tmp2;
	spinor plus;

	hmc_complex twistfactor = {1., mubar_in};
	hmc_complex twistfactor_minus = {1., -1.*mubar_in};

	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size) {

		/** @todo this must be done more efficient */
		st_index pos = (id_tmp % 2 == 0) ? get_even_site(id_tmp / 2) : get_odd_site(id_tmp / 2);

		out_tmp = set_spinor_zero();
		out_tmp2 = set_spinor_zero();
		//get input spinor
		plus = get_spinor_from_field(in, pos.space, pos.time);
		//Diagonalpart: this is just the normal tm-diagonal matrix
		out_tmp = M_diag_tm_local(plus, twistfactor, twistfactor_minus);
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
