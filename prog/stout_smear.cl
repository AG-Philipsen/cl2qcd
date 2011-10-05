/**
 * @file stout-smearing of the gaugefield
 */

//this is done after the tmlqcd-analogue (stout_smear.c)
__kernel void stout_smear(int   iterations, hmc_float rho, __global ocl_s_gaugefield * u_inout)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	int n, t, mu;

	/** @todo no saving of the intermediate fields implemented!! This should be necessary for the force-calculation*/
	/** @todo how to get rho_iter and mu in here? As arguments? */
	for(int id_tmp = id; id_tmp < VOL4D; id_tmp += global_size) {
		/** @todo this must be done more efficient */
		if(id_tmp < VOLSPACE * NTIME / 2)
			get_even_site(id_tmp, &n, &t);
		else
			get_odd_site(id_tmp - (VOLSPACE * NTIME / 2), &n, &t);

		//calc staple
		Matrix3x3 wk_staple = calc_staple(u_inout, n, t, mu) ;
		//staple*rho
		wk_staple = multiply_matrix3x3_by_real(wk_staple, rho);

		//omega = staple * u^dagger
		Matrixsu3 gauge_local;
		gauge_local = get_matrixsu3(u_inout, n, t, mu);
		Matrix3x3 tmp2 = matrix_su3to3x3(gauge_local);
		Matrix3x3 omega = multiply_matrix3x3_dagger(wk_staple, tmp2);

		//project out anti-hermitian traceless part
		//reuse gauge_local here again
		gauge_local = project_anti_herm(omega);

		//exponentiate
		ae p;
		p = tr_lambda_u(omega) ;
		// -2.0 to get su3 to su3adjoint consistency
		p = ae_times_factor(p, -0.5);

		Matrixsu3 Exp_p = build_su3matrix_by_exponentiation((p), 1);

		//new_gauge_local = Exp_p * gauge_local */
		Matrixsu3 new_gauge_local = multiply_matrixsu3(Exp_p, gauge_local);
		put_matrixsu3(u_inout, new_gauge_local, n, t, mu);

	}
}
