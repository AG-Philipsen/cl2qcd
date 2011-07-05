/** @file
 * Device code for the heatbath overrelaxation
 */

void inline perform_overrelaxing(__global ocl_s_gaugefield* gaugefield, const hmc_float beta, const int mu, __global hmc_ocl_ran * rnd, int pos, int t, int id)
{

	Matrixsu3 U;
	Matrix3x3 W;
	Matrix3x3 staplematrix;

	Matrixsu2 w;
	Matrixsu2_pauli w_pauli;
	hmc_float k;
	int order[3];

	random_1_2_3(order, &rnd[id]);
	U = get_matrixsu3(gaugefield, pos, t, mu);
	U = project_su3(U);

	staplematrix = calc_staple(gaugefield, pos, t, mu);

	Matrixsu3 extW;

	for(int i = 0; i < NC; i++) {
		W = matrix_su3to3x3 (U);
		W = multiply_matrix3x3 (W, staplematrix);

		w = reduction(W, order[i]);

		w_pauli.e00 = 0.5 * (w.e00.re + w.e11.re);
		w_pauli.e01 = 0.5 * (w.e01.im + w.e10.im);
		w_pauli.e10 = 0.5 * (w.e01.re - w.e10.re);
		w_pauli.e11 = 0.5 * (w.e00.im - w.e11.im);
		k = sqrt(  w_pauli.e00 * w_pauli.e00 +  w_pauli.e01 * w_pauli.e01 + w_pauli.e10 * w_pauli.e10 + w_pauli.e11 * w_pauli.e11  );

		w.e00.re = (w_pauli.e00 * w_pauli.e00 - w_pauli.e01 * w_pauli.e01 - w_pauli.e10 * w_pauli.e10 - w_pauli.e11 * w_pauli.e11) / k / k;
		w.e00.im = (-2.*w_pauli.e00 * w_pauli.e11) / k / k;
		w.e01.re = (-2.*w_pauli.e00 * w_pauli.e10) / k / k;
		w.e01.im = (-2.*w_pauli.e00 * w_pauli.e01) / k / k;
		w.e10.re = (2.*w_pauli.e00 * w_pauli.e10) / k / k;
		w.e10.im = (-2.*w_pauli.e00 * w_pauli.e01) / k / k;
		w.e11.re = (w_pauli.e00 * w_pauli.e00 - w_pauli.e01 * w_pauli.e01 - w_pauli.e10 * w_pauli.e10 - w_pauli.e11 * w_pauli.e11) / k / k;
		w.e11.im = (2.*w_pauli.e00 * w_pauli.e11) / k / k;

		extW = extend (order[i], w);
		U = multiply_matrixsu3(extW, U);
	}

	put_matrixsu3(gaugefield, U, pos, t, mu);
}

__kernel void overrelax_even(__global ocl_s_gaugefield* gaugefield, const hmc_float beta, const int mu, __global hmc_ocl_ran * rnd)
{
	int t, pos, id, id_tmp, size;
	id_tmp = get_global_id(0);
	size = get_global_size(0);
	for(id = id_tmp; id < VOLSPACE * NTIME / 2; id += size) {
		get_even_site(id, &pos, &t);
		perform_overrelaxing(gaugefield, beta, mu, rnd, pos, t, id_tmp);
	}
}

__kernel void overrelax_odd(__global ocl_s_gaugefield* gaugefield, const hmc_float beta, const int mu, __global hmc_ocl_ran * rnd)
{
	int t, pos, id, id_tmp, size;
	id_tmp = get_global_id(0);
	size = get_global_size(0);
	for(id = id_tmp; id < VOLSPACE * NTIME / 2; id += size) {
		get_odd_site(id, &pos, &t);
		perform_overrelaxing(gaugefield, beta, mu, rnd, pos, t, id_tmp);
	}
}


