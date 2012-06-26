inline void gauge_force_tlsym_per_link(__global const Matrixsu3StorageType * const restrict field, __global aeStorageType * const restrict out, const st_index pos, const dir_idx dir) {
		Matrix3x3 V = calc_rectangles_staple(field, pos.space, pos.time, dir);
		Matrixsu3 U = get_matrixsu3(field, pos.space, pos.time, dir);
		V = multiply_matrix3x3 (matrix_su3to3x3(U), V);
		ae out_tmp = tr_lambda_u(V);

		hmc_float factor = -C1 * BETA / 3.;
		int global_link_pos = get_global_link_pos(dir, pos.space, pos.time);
		update_gaugemomentum(out_tmp, factor , global_link_pos, out);
}

#ifdef _USEGPU_
__attribute__((reqd_work_group_size(128, 1, 1)))
#endif
__kernel void gauge_force_tlsym(__global const Matrixsu3StorageType * const restrict field, __global aeStorageType * const restrict out)
{
#ifndef _USE_RECT_
	//this kernel should not be called if rectangles are not activated
#error should not be used if rect is deactivated
#endif

	//tlSym improved Gauge force is factor*Im(i Tr(T_i U V))
	//   with T_i being the SU3-Generator in i-th direction and V the staplematrix
	//   and the factor being 0 (for standard Wilson-action) and -c1 * beta / NC (for tlSym)
	PARALLEL_FOR(id_tmp, VOL4D / 2 * NDIM) {
		//calc link-pos and mu out of the index
		//NOTE: this is not necessarily equal to the geometric  conventions, one just needs a one-to-one correspondence between thread-id and (n,t,mu) here
		int2 pos_tmp;
		pos_tmp.x = id_tmp % (VOL4D / 2);
		pos_tmp.y = id_tmp / (VOL4D / 2);
		st_index pos = get_even_site(pos_tmp.x);
		gauge_force_tlsym_per_link(field, out, pos, pos_tmp.y);
		pos = get_odd_site(pos_tmp.x);
		gauge_force_tlsym_per_link(field, out, pos, pos_tmp.y);
	}
}
