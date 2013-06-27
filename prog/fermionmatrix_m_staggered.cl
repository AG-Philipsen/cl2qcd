/** @file
 * M normal staggered fermionmatrix
 *
 * In order to clarify notation and not to confuse things, let us recall some expressions.
 * In the staggered formulation, after the "staggering process", we remain with one component
 * field per site (an su3vec per site if we take into account the colour degree of freedom). 
 * Since in an HMC algorithm we do not deal directly with such one-component field (in fact,
 * first we integrate over the fermionic field getting the fermionic determinant and then
 * we introduce a pseudofermionic field to rewrite the fermionic determinat in a suitable way),
 * here we will focus solely on the Dirac operator that we will indicate with M. Referring for
 * instance to the Rothe (Lattice Gauge Theories, Eq.(6.26), pag 94), we have:
 *
 * \f[
   M_{n,m}=\frac{1}{2} \sum_{\mu=0}^3 \eta_\mu(n)\Bigl[U_\mu(n)\,\delta_{n+\hat\mu,m}
                                                     - U^\dag_\mu(n-\hat\mu)\,\delta_{n-\hat\mu,m}\Bigr]
                                                     + m_0\,\delta_{n,m}
   \f]
 *
 * \internal   
 *
 *               ___            _                                                     _ 
 *            1  \             |                                                       |   
 * M_{n,m} = --- /__  eta_mu(n)| U_mu(n) delta_{n+mu,m} - U^dag_mu(n-mu) delta_{n-mu,m}| + m delta_{n,m} 
 *            2   mu           |_                                                     _|   
 *
 *
 * We refer to M as M=m+D_KS
 *
 * \endinternal
 * 
 * Here in this kernel we make the M operator act on the field "in" getting the result in the
 * field "out". Note that the field "field" is the gauge field (i.e. the link variables) that
 * are needed in the functions dslash_local_x.
 *
 * @note This kernel takes into account the staggered phases eta_mu(n), since they
 *       will be NOT included in links (as in the phi-algorithm from Gottlieb and Toussaint).
 *       The reason why we make such a choice is that, on GPU, it should be more efficient to
 *       calculate staggered phases whenever we need them, instead of including them in
 *       links and hence projecting links to SU(3) group sometimes (one should do such
 *       a projection because if eta_mu is -1 then you go out from SU(3) group).
 *       Observe that, avoiding such projection, we get rid of this approximation. Actually
 *       this is not completely true because there is still a projection to be done because
 *       of the numerical errors in the integration of the molecular dinamics equations.
 * \par
 * @note In this kernel \em no \em even-odd \em preconditioning is assumed: vectors "in" and "out" have
 *       VOL4D components.
 */


__kernel void M_staggered(__global const su3vec * const restrict in, __global const Matrixsu3StorageType * const restrict field, __global su3vec * const restrict out, hmc_float mass_in)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	su3vec out_tmp, out_tmp2;

	for(int id_local = id; id_local < SPINORFIELDSIZE_LOCAL; id_local += global_size) {

		/** @todo this must be done more efficient */
		st_index pos = (id_local % 2 == 0) ? get_even_st_idx_local(id_local / 2) : get_odd_st_idx_local(id_local / 2);

		
		//From now on we adopt the notation M = m + D_KS
		
		//Diagonal part: m * in(n)
		out_tmp = get_su3vec_from_field(in, pos.space, pos.time);
		out_tmp = su3vec_times_real(out_tmp, mass_in);
		
		//Non-diagonal part: calc D_KS
		out_tmp2 = D_KS_local(in, field, pos, TDIR);
		out_tmp = su3vec_acc(out_tmp, out_tmp2);
		out_tmp2 = D_KS_local(in, field, pos, XDIR);
		out_tmp = su3vec_acc(out_tmp, out_tmp2);
		out_tmp2 = D_KS_local(in, field, pos, YDIR);
		out_tmp = su3vec_acc(out_tmp, out_tmp2);
		out_tmp2 = D_KS_local(in, field, pos, ZDIR);
		out_tmp = su3vec_acc(out_tmp, out_tmp2);
		
		put_su3vec_to_field(out_tmp, out, pos.space, pos.time);
	}
}

#if 0
 // The following kernel uses the d_slash_local_x kernels and it was the first version written.
 // At the moment we leave it out. Uncomment it out if needed for some reason.

__kernel void M_staggered_old(__global const su3vec * const restrict in, __global const Matrixsu3StorageType * const restrict field, __global su3vec * const restrict out, hmc_float mass_in)
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

		
		//From now on we adopt the notation M = m + D_KS
		//Diagonal part: m * in(n)
		out_tmp = get_su3vec_from_field(in, pos.space, pos.time);
		out_tmp = su3vec_times_real(out_tmp, mass_in);
		
		//Non-diagonal part: calc D_KS
		out_tmp2 = dslash_local_0(in, field, pos.space, pos.time);
		out_tmp = su3vec_acc(out_tmp, out_tmp2);
		out_tmp2 = dslash_local_1(in, field, pos.space, pos.time);
		out_tmp = su3vec_acc(out_tmp, out_tmp2);
		out_tmp2 = dslash_local_2(in, field, pos.space, pos.time);
		out_tmp = su3vec_acc(out_tmp, out_tmp2);
		out_tmp2 = dslash_local_3(in, field, pos.space, pos.time);
		out_tmp = su3vec_acc(out_tmp, out_tmp2);
				
		put_su3vec_to_field(out_tmp, out, pos.space, pos.time);
	}
}
#endif
