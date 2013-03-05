inline void gauge_force_tlsym_per_link(__global const Matrixsu3StorageType * const restrict field, __global aeStorageType * const restrict out, const st_index pos, const dir_idx dir) {
		Matrix3x3 V = calc_rectangles_staple(field, pos.space, pos.time, dir);
		Matrixsu3 U = get_matrixsu3(field, pos.space, pos.time, dir);
		V = multiply_matrix3x3 (matrix_su3to3x3(U), V);
		ae out_tmp = tr_lambda_u(V);

		hmc_float factor = -C1 * BETA / 3.;
		int global_link_pos = get_link_idx(dir, pos);
		update_gaugemomentum(out_tmp, factor , global_link_pos, out);
}

__kernel void gauge_force_tlsym(__global const Matrixsu3StorageType * const restrict field, __global aeStorageType * const restrict out)
{
#ifndef _USE_RECT_
	//this kernel should not be called if rectangles are not activated
#error should not be used if rect is deactivated
#endif

	//tlSym improved Gauge force is factor*Im(i Tr(T_i U V))
	//   with T_i being the SU3-Generator in i-th direction and V the staplematrix
	//   and the factor being 0 (for standard Wilson-action) and -c1 * beta / NC (for tlSym)
	PARALLEL_FOR(id_tmp, VOL4D * NDIM) {
		//calc link-pos and mu out of the index
		//NOTE: this is not necessarily equal to the geometric  conventions, one just needs a one-to-one correspondence between thread-id and (n,t,mu) here
		const size_t pos_tmp = id_tmp % VOL4D;
		const size_t dir     = id_tmp / VOL4D;
		const st_index pos = (pos_tmp >= VOL4D / 2) ? get_even_site(pos_tmp - (VOL4D / 2)) : get_odd_site(pos_tmp);
		gauge_force_tlsym_per_link(field, out, pos, dir);
	}
}

/*
 A completely different variant
*/

__kernel void gauge_force_tlsym_multipass1_tpe(__global const Matrixsu3StorageType * const restrict field, __global Matrix3x3StorageType * const restrict tmp)
{
	const size_t id_tmp = get_global_id(0);
	if(id_tmp < VOL4D * NDIM) {
		//calc link-pos and mu out of the index
		//NOTE: this is not necessarily equal to the geometric  conventions, one just needs a one-to-one correspondence between thread-id and (n,t,mu) here
		const size_t pos_tmp = id_tmp % VOL4D;
		const size_t dir     = id_tmp / VOL4D;
		const st_index pos = (pos_tmp >= VOL4D / 2) ? get_even_site(pos_tmp - (VOL4D / 2)) : get_odd_site(pos_tmp);

		Matrix3x3 staple = zero_matrix3x3();
#ifdef _USEGPU_
#pragma unroll 3 // unroll required for proper register reuse when using newer Catalysts on Cypress
#endif
		for(int i = 0; i < NDIM - 1 ; i++) {
			int nu = (dir + i + 1) % NDIM;
			staple = add_matrix3x3(staple, local_rectangles_staple_1(field, pos.space, pos.time, dir, nu));
		}
		put_matrix3x3(tmp, staple, pos.space, pos.time, dir);
	}
}

__kernel void gauge_force_tlsym_multipass2_tpe(__global const Matrixsu3StorageType * const restrict field, __global Matrix3x3StorageType * const restrict tmp)
{
	const size_t id_tmp = get_global_id(0);
	if(id_tmp < VOL4D * NDIM) {
		//calc link-pos and mu out of the index
		//NOTE: this is not necessarily equal to the geometric  conventions, one just needs a one-to-one correspondence between thread-id and (n,t,mu) here
		const size_t pos_tmp = id_tmp % VOL4D;
		const size_t dir     = id_tmp / VOL4D;
		const st_index pos = (pos_tmp >= VOL4D / 2) ? get_even_site(pos_tmp - (VOL4D / 2)) : get_odd_site(pos_tmp);

		Matrix3x3 staple = get_matrix3x3(tmp, pos.space, pos.time, dir);
#ifdef _USEGPU_
#pragma unroll 3 // unroll required for proper register reuse when using newer Catalysts on Cypress
#endif
		for(int i = 0; i < NDIM - 1 ; i++) {
			int nu = (dir + i + 1) % NDIM;
			staple = add_matrix3x3(staple, local_rectangles_staple_2(field, pos.space, pos.time, dir, nu));
		}
		put_matrix3x3(tmp, staple, pos.space, pos.time, dir);
	}
}

__kernel void gauge_force_tlsym_multipass3_tpe(__global const Matrixsu3StorageType * const restrict field, __global Matrix3x3StorageType * const restrict tmp)
{
	const size_t id_tmp = get_global_id(0);
	if(id_tmp < VOL4D * NDIM) {
		//calc link-pos and mu out of the index
		//NOTE: this is not necessarily equal to the geometric  conventions, one just needs a one-to-one correspondence between thread-id and (n,t,mu) here
		const size_t pos_tmp = id_tmp % VOL4D;
		const size_t dir     = id_tmp / VOL4D;
		const st_index pos = (pos_tmp >= VOL4D / 2) ? get_even_site(pos_tmp - (VOL4D / 2)) : get_odd_site(pos_tmp);

		Matrix3x3 staple = get_matrix3x3(tmp, pos.space, pos.time, dir);
#ifdef _USEGPU_
#pragma unroll 3 // unroll required for proper register reuse when using newer Catalysts on Cypress
#endif
		for(int i = 0; i < NDIM - 1 ; i++) {
			int nu = (dir + i + 1) % NDIM;
			staple = add_matrix3x3(staple, local_rectangles_staple_3(field, pos.space, pos.time, dir, nu));
		}
		put_matrix3x3(tmp, staple, pos.space, pos.time, dir);
	}
}

__kernel void gauge_force_tlsym_multipass4_tpe(__global const Matrixsu3StorageType * const restrict field, __global Matrix3x3StorageType * const restrict tmp)
{
	const size_t id_tmp = get_global_id(0);
	if(id_tmp < VOL4D * NDIM) {
		//calc link-pos and mu out of the index
		//NOTE: this is not necessarily equal to the geometric  conventions, one just needs a one-to-one correspondence between thread-id and (n,t,mu) here
		const size_t pos_tmp = id_tmp % VOL4D;
		const size_t dir     = id_tmp / VOL4D;
		const st_index pos = (pos_tmp >= VOL4D / 2) ? get_even_site(pos_tmp - (VOL4D / 2)) : get_odd_site(pos_tmp);

		Matrix3x3 staple = get_matrix3x3(tmp, pos.space, pos.time, dir);
#ifdef _USEGPU_
#pragma unroll 3 // unroll required for proper register reuse when using newer Catalysts on Cypress
#endif
		for(int i = 0; i < NDIM - 1 ; i++) {
			int nu = (dir + i + 1) % NDIM;
			staple = add_matrix3x3(staple, local_rectangles_staple_4(field, pos.space, pos.time, dir, nu));
		}
		put_matrix3x3(tmp, staple, pos.space, pos.time, dir);
	}
}

__kernel void gauge_force_tlsym_multipass5_tpe(__global const Matrixsu3StorageType * const restrict field, __global Matrix3x3StorageType * const restrict tmp)
{
	const size_t id_tmp = get_global_id(0);
	if(id_tmp < VOL4D * NDIM) {
		//calc link-pos and mu out of the index
		//NOTE: this is not necessarily equal to the geometric  conventions, one just needs a one-to-one correspondence between thread-id and (n,t,mu) here
		const size_t pos_tmp = id_tmp % VOL4D;
		const size_t dir     = id_tmp / VOL4D;
		const st_index pos = (pos_tmp >= VOL4D / 2) ? get_even_site(pos_tmp - (VOL4D / 2)) : get_odd_site(pos_tmp);

		Matrix3x3 staple = get_matrix3x3(tmp, pos.space, pos.time, dir);
#ifdef _USEGPU_
#pragma unroll 3 // unroll required for proper register reuse when using newer Catalysts on Cypress
#endif
		for(int i = 0; i < NDIM - 1 ; i++) {
			int nu = (dir + i + 1) % NDIM;
			staple = add_matrix3x3(staple, local_rectangles_staple_5(field, pos.space, pos.time, dir, nu));
		}
		put_matrix3x3(tmp, staple, pos.space, pos.time, dir);
	}
}

__kernel void gauge_force_tlsym_multipass6_tpe(__global const Matrixsu3StorageType * const restrict field, __global aeStorageType * const restrict out, __global Matrix3x3StorageType * const restrict tmp)
{
	const size_t id_tmp = get_global_id(0);
	if(id_tmp < VOL4D * NDIM) {
		//calc link-pos and mu out of the index
		//NOTE: this is not necessarily equal to the geometric  conventions, one just needs a one-to-one correspondence between thread-id and (n,t,mu) here
		const size_t pos_tmp = id_tmp % VOL4D;
		const size_t dir     = id_tmp / VOL4D;
		const st_index pos = (pos_tmp >= VOL4D / 2) ? get_even_site(pos_tmp - (VOL4D / 2)) : get_odd_site(pos_tmp);

		Matrix3x3 staple = get_matrix3x3(tmp, pos.space, pos.time, dir);
#ifdef _USEGPU_
#pragma unroll 3 // unroll required for proper register reuse when using newer Catalysts on Cypress
#endif
		for(int i = 0; i < NDIM - 1 ; i++) {
			int nu = (dir + i + 1) % NDIM;
			staple = add_matrix3x3(staple, local_rectangles_staple_6(field, pos.space, pos.time, dir, nu));
		}

		Matrixsu3 U = get_matrixsu3(field, pos.space, pos.time, dir);
		staple = multiply_matrix3x3 (matrix_su3to3x3(U), staple);
		ae out_tmp = tr_lambda_u(staple);

		hmc_float factor = -C1 * BETA / 3.;
		int global_link_pos = get_link_idx(dir, pos);
		update_gaugemomentum(out_tmp, factor , global_link_pos, out);
	}
}
