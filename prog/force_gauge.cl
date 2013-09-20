/** @file
 * Kernel for the standard Wilson-action gauge force calculation.
 * 
 * According to the Gattringer book (page 198) the gauge part of the force should be
 * @code
 * F_G = - \beta/6 * \sum_{k=1}^8 T_k Tr[i * T_k * (U*V - V^\dag*U^\dag)] 
 * 
 *     = - \beta/12 * i * (U*V - V^\dag*U^\dag)
 * @endcode
 * where V are the well known staples. Nevertheless this expression is not in
 * agreement with the Gottlieb and Toussaint work (Hybrid-molecular dynamics
 * algorithm for numerical simulation of QCD), where we found
 * @code
 * F_G =  - \beta/6 * i * (U*V - V^\dag*U^\dag) .
 * @endcode
 * Then there is a factor 1/2 of discrepancy.
 * 
 * @note It is worth recalling that an overall factor in the force calculation does
 *       not affect the correctness of the (R)HMC algorithm. In fact (thinking to
 *       leapfrog, see page 197 Gattringer), this factor would mean to have a different
 *       time step size in the integration of the equations of motion. Of course if one
 *       wants to make the fields evolve for a fixed time in, let's say, 100 steps,
 *       he will deduce a certain time step. A missing overall factor in the force
 *       calculation would affect the time step and therefore all the fields would
 *       evolve for less or more time!
 * 
 * @attention Now, here, as one can see reading the code, the matrix U*V is calculated
 *            and then the function tr_lambda_u is called on it. This function returns
 *            Tr[i * T_k * (U*V - V^\dag*U^\dag)] that is calculated ONLY by @f$ \beta/3 @f$.
 *            We are then coeherent with the Gottlieb and Toussaint result.
 */

inline void gauge_force_per_link(__global const Matrixsu3StorageType * const restrict field, __global aeStorageType * const restrict out, const st_index pos, const dir_idx dir) {
	Matrix3x3 V = calc_staple(field, pos.space, pos.time, dir);
	Matrixsu3 U = get_matrixsu3(field, pos.space, pos.time, dir);
	V = multiply_matrix3x3 (matrix_su3to3x3(U), V);
	ae out_tmp = tr_lambda_u(V);

	hmc_float factor = -BETA / 3.;
#ifdef _USE_RECT_
	factor = factor * C0;
#endif
	int global_link_pos = get_link_idx(dir, pos);
	update_gaugemomentum(out_tmp, factor , global_link_pos, out);
}

__kernel void gauge_force(__global const Matrixsu3StorageType * const restrict field, __global aeStorageType * const restrict out)
{
	//Gauge force is factor*Im(i Tr(T_i U V))
	//   with T_i being the SU3-Generator in i-th direction and V the staplematrix
	//   and the factor being - beta / NC (for standard Wilson-action) and -c0 * beta / NC (for tlSym)
	PARALLEL_FOR(id_local, VOL4D_LOCAL * NDIM) {
		//calc link-pos and mu out of the index
		//NOTE: this is not necessarily equal to the geometric  conventions, one just needs a one-to-one correspondence between thread-id and (n,t,mu) here
		const int pos_local = id_local % VOL4D_LOCAL;
		const int dir     = id_local / VOL4D_LOCAL;
		const st_index pos = (pos_local >= VOL4D_LOCAL / 2) ? get_even_st_idx_local(pos_local - (VOL4D_LOCAL / 2)) : get_odd_st_idx_local(pos_local);

		gauge_force_per_link(field, out, pos, dir);
	}
}
