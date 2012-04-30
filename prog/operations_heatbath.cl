/** @file
 * Device code for the heatbath update
 */

//operations_heatbath.cl

Matrixsu2_pauli SU2Update(const hmc_float alpha, prng_state * const restrict rnd)
{
	Matrixsu2_pauli out;

	hmc_float delta;
	hmc_float a0 ;
	hmc_float eta ;
	do {
		double4 rands = prng_double4(rnd);
		delta = -log(rands.x) / alpha * pow(cos((hmc_float)(2.f * PI * rands.y)), (hmc_float) 2.f) - log(rands.z) / alpha;
		a0 = 1. - delta;
		eta = rands.w;
	} while ( (1. - 0.5 * delta) < eta * eta);
	prng_synchronize(rnd);

	double4 rands = prng_double4(rnd);
	hmc_float phi = 2.*PI * rands.x;
	hmc_float theta = asin((hmc_float)(2.f * rands.y - 1.f));
	out.e00 = a0;
	out.e01 = sqrt(1.f - a0 * a0) * cos(theta) * cos(phi);
	out.e10 = sqrt(1.f - a0 * a0) * cos(theta) * sin(phi);
	out.e11 = sqrt(1.f - a0 * a0) * sin(theta);

	return out;
}

void inline perform_heatbath(__global Matrixsu3StorageType * const restrict dest_gaugefield, __global const Matrixsu3StorageType * const restrict src_gaugefield, const int mu, prng_state * const restrict rnd, const int pos, const int t)
{
	Matrix3x3 staplematrix;

#ifdef _ANISO_
	//Compute staple, comprises whole anisotropy
	if (mu == 0) {
		staplematrix = calc_staple(src_gaugefield, pos, t, mu);
		staplematrix = multiply_matrix3x3_by_real (staplematrix, XI_0 );
	}

	else {
		Matrix3x3 staplematrix_sigma;
		Matrix3x3 staplematrix_tau;
		staplematrix_sigma = calc_staple_sigma(src_gaugefield, pos, t, mu);
		staplematrix_sigma = multiply_matrix3x3_by_real (staplematrix_sigma, 1 / XI_0 );
		staplematrix_tau = calc_staple_tau(src_gaugefield, pos, t, mu);
		staplematrix_tau = multiply_matrix3x3_by_real (staplematrix_tau, XI_0 );
		staplematrix = add_matrix3x3 ( staplematrix_sigma, staplematrix_tau );
	}
#else
	staplematrix = calc_staple(src_gaugefield, pos, t, mu);
#endif

	Matrixsu3 U = get_matrixsu3(src_gaugefield, pos, t, mu);
	U = project_su3(U);

	int3 order = prng_123(rnd);

	for(int i = 0; i < NC; i++) {
		int order_i = (i == 0) ? order.x : ((i == 1) ? order.y : order.z);

		Matrix3x3 W = matrix_su3to3x3 (U);
		W = multiply_matrix3x3 (W, staplematrix);

		Matrixsu2 w = reduction(W, order_i);

		Matrixsu2_pauli w_pauli;
		w_pauli.e00 = 0.5 * (w.e00.re + w.e11.re);
		w_pauli.e01 = 0.5 * (w.e01.im + w.e10.im);
		w_pauli.e10 = 0.5 * (w.e01.re - w.e10.re);
		w_pauli.e11 = 0.5 * (w.e00.im - w.e11.im);
		hmc_float k = sqrt(  w_pauli.e00 * w_pauli.e00 +  w_pauli.e01 * w_pauli.e01 + w_pauli.e10 * w_pauli.e10 + w_pauli.e11 * w_pauli.e11  );

		hmc_float beta_new =  2.*BETA / NC * k;

		Matrixsu2_pauli r_pauli = SU2Update(beta_new, rnd);

		w.e00.re = (r_pauli.e00 * w_pauli.e00 + r_pauli.e01 * w_pauli.e01 + r_pauli.e10 * w_pauli.e10 + r_pauli.e11 * w_pauli.e11 ) / k;
		w.e00.im = (w_pauli.e00 * r_pauli.e11 - w_pauli.e11 * r_pauli.e00 + r_pauli.e01 * w_pauli.e10 - r_pauli.e10 * w_pauli.e01 ) / k;
		w.e01.re = (w_pauli.e00 * r_pauli.e10 - w_pauli.e10 * r_pauli.e00 + r_pauli.e11 * w_pauli.e01 - r_pauli.e01 * w_pauli.e11 ) / k;
		w.e01.im = (w_pauli.e00 * r_pauli.e01 - w_pauli.e01 * r_pauli.e00 + r_pauli.e10 * w_pauli.e11 - r_pauli.e11 * w_pauli.e10 ) / k;
		w.e10.re = -(w_pauli.e00 * r_pauli.e10 - w_pauli.e10 * r_pauli.e00 + r_pauli.e11 * w_pauli.e01 - r_pauli.e01 * w_pauli.e11 ) / k;
		w.e10.im = (w_pauli.e00 * r_pauli.e01 - w_pauli.e01 * r_pauli.e00 + r_pauli.e10 * w_pauli.e11 - r_pauli.e11 * w_pauli.e10 ) / k;
		w.e11.re = (r_pauli.e00 * w_pauli.e00 + r_pauli.e01 * w_pauli.e01 + r_pauli.e10 * w_pauli.e10 + r_pauli.e11 * w_pauli.e11 ) / k;
		w.e11.im = -(w_pauli.e00 * r_pauli.e11 - w_pauli.e11 * r_pauli.e00 + r_pauli.e01 * w_pauli.e10 - r_pauli.e10 * w_pauli.e01 ) / k;

		Matrixsu3 extW = extend (order_i, w);


		U = multiply_matrixsu3 (extW, U);
	}

	put_matrixsu3(dest_gaugefield, U, pos, t, mu);

}

void inline perform_overrelaxing(__global Matrixsu3StorageType * const restrict dest_gaugefield, __global const Matrixsu3StorageType * const restrict src_gaugefield, const int mu, prng_state * const restrict rnd, const int pos, const int t)
{
	Matrix3x3 staplematrix;
#ifdef _ANISO_
	//Compute staple, comprises whole anisotropy
	if (mu == 0) {
		staplematrix = calc_staple(src_gaugefield, pos, t, mu);
		staplematrix = multiply_matrix3x3_by_real (staplematrix, XI_0 );
	}

	else {
		Matrix3x3 staplematrix_sigma;
		Matrix3x3 staplematrix_tau;
		staplematrix_sigma = calc_staple_sigma(src_gaugefield, pos, t, mu);
		staplematrix_sigma = multiply_matrix3x3_by_real (staplematrix_sigma, 1 / XI_0 );
		staplematrix_tau = calc_staple_tau(src_gaugefield, pos, t, mu);
		staplematrix_tau = multiply_matrix3x3_by_real (staplematrix_tau, XI_0 );
		staplematrix = add_matrix3x3 ( staplematrix_sigma, staplematrix_tau );
	}
#else
	staplematrix = calc_staple(src_gaugefield, pos, t, mu);
#endif
	Matrixsu3 U = get_matrixsu3(src_gaugefield, pos, t, mu);
	U = project_su3(U);

	int3 order = prng_123(rnd);

	for(int i = 0; i < NC; i++) {
		int order_i = (i == 0) ? order.x : ((i == 1) ? order.y : order.z);
		Matrix3x3 W = matrix_su3to3x3 (U);
		W = multiply_matrix3x3 (W, staplematrix);

		Matrixsu2 w = reduction(W, order_i);

		Matrixsu2_pauli w_pauli;
		w_pauli.e00 = 0.5 * (w.e00.re + w.e11.re);
		w_pauli.e01 = 0.5 * (w.e01.im + w.e10.im);
		w_pauli.e10 = 0.5 * (w.e01.re - w.e10.re);
		w_pauli.e11 = 0.5 * (w.e00.im - w.e11.im);
		hmc_float k = sqrt(  w_pauli.e00 * w_pauli.e00 +  w_pauli.e01 * w_pauli.e01 + w_pauli.e10 * w_pauli.e10 + w_pauli.e11 * w_pauli.e11  );

		w.e00.re = (w_pauli.e00 * w_pauli.e00 - w_pauli.e01 * w_pauli.e01 - w_pauli.e10 * w_pauli.e10 - w_pauli.e11 * w_pauli.e11) / k / k;
		w.e00.im = (-2.*w_pauli.e00 * w_pauli.e11) / k / k;
		w.e01.re = (-2.*w_pauli.e00 * w_pauli.e10) / k / k;
		w.e01.im = (-2.*w_pauli.e00 * w_pauli.e01) / k / k;
		w.e10.re = (2.*w_pauli.e00 * w_pauli.e10) / k / k;
		w.e10.im = (-2.*w_pauli.e00 * w_pauli.e01) / k / k;
		w.e11.re = (w_pauli.e00 * w_pauli.e00 - w_pauli.e01 * w_pauli.e01 - w_pauli.e10 * w_pauli.e10 - w_pauli.e11 * w_pauli.e11) / k / k;
		w.e11.im = (2.*w_pauli.e00 * w_pauli.e11) / k / k;

		Matrixsu3 extW = extend (order_i, w);
		U = multiply_matrixsu3(extW, U);
	}

	put_matrixsu3(dest_gaugefield, U, pos, t, mu);
}
