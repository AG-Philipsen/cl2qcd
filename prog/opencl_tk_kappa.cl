/** @file
 * Device code for the kappa measurement
 */

//opencl_tk_kappa.cl

__kernel void kappa_karsch_gpu(__global Matrixsu3StorageType * gaugefield, const hmc_float beta, __global hmc_float * kappa_karsch_val)
{

	//Compute diagonal spatial components of the energy-momentum-tensor
	hmc_float tdiag_11 [VOL4D];
	hmc_float tdiag_22 [VOL4D];
	hmc_float tdiag_33 [VOL4D];
	//a = 2T_11 - T_22 _ T_33
	hmc_float a [VOL4D];
	//b = 2T_22 - T_11 _ T_33
	hmc_float b [VOL4D];
	//c = 2T_33 - T_22 _ T_11
	hmc_float c [VOL4D];

	for (int t = 0; t < NTIME; t++) {
		for (int n = 0; n < VOLSPACE; n++) {
			//Compute required plaquettes
			Matrixsu3 temp;

			temp = local_plaquette(gaugefield, n, t, 1, 0);
			hmc_float plaq_10 = trace_matrixsu3(temp).re;
			temp = local_plaquette(gaugefield, n, t, 2, 0);
			hmc_float plaq_20 = trace_matrixsu3(temp).re;
			temp = local_plaquette(gaugefield, n, t, 3, 0);
			hmc_float plaq_30 = trace_matrixsu3(temp).re;
			temp = local_plaquette(gaugefield, n, t, 1, 2);
			hmc_float plaq_12 = trace_matrixsu3(temp).re;
			temp = local_plaquette(gaugefield, n, t, 1, 3);
			hmc_float plaq_13 = trace_matrixsu3(temp).re;
			temp = local_plaquette(gaugefield, n, t, 3, 2);
			hmc_float plaq_32 = trace_matrixsu3(temp).re;

			int point = n + VOLSPACE * t;

			tdiag_11 [point] = plaq_10 + plaq_12 + plaq_13 - plaq_20 - plaq_30 - plaq_32;
			tdiag_22 [point] = plaq_20 + plaq_12 + plaq_32 - plaq_10 - plaq_30 - plaq_13;
			tdiag_33 [point] = plaq_30 + plaq_13 + plaq_32 - plaq_10 - plaq_20 - plaq_12;

			a[point] = 2.0 * tdiag_11[point] - tdiag_22[point] - tdiag_33[point];
			b[point] = 2.0 * tdiag_22[point] - tdiag_11[point] - tdiag_33[point];
			c[point] = 2.0 * tdiag_33[point] - tdiag_22[point] - tdiag_11[point];
		}
	}

	hmc_float deltak = 2 * PI / (hmc_float) NSPACE;
	hmc_float result = 0.0;


	for (int x_3 = 0; x_3 < NSPACE; x_3++) {
		for (int y_3 = 0; y_3 < x_3; y_3++) {
			hmc_float factor = 1.0 - cos(deltak * (hmc_float) (x_3 - y_3));
			for (int x_t = 0; x_t < NTIME; x_t++) {
				for (int y_t = 0; y_t < NTIME; y_t++) {
					uint3 coord_y, coord_y;
					for (int coord_x.x = 0; coord_x.x < NSPACE; coord_x.x++) {
						for (int coord_y.x = 0; coord_y.x < NSPACE; coord_y.x++) {
							for (int coord_x.y = 0; coord_x.y < NSPACE; coord_x.y++) {
								for (int coord_y.y = 0; coord_y.y < NSPACE; coord_y.y++) {

									//new method get_tnspace which gives n_x+VOLSPACE*x_t
									int n_x = get_nspace (coord_x);
									int point_x = n_x + VOLSPACE * x_t;
									int n_y = get_nspace (coord_y);
									int point_y = n_y + VOLSPACE * y_t;

									result += factor * ( tdiag_11 [point_x] * a[point_y]
									                     + tdiag_22 [point_x] * b[point_y]
									                     + tdiag_33 [point_x] * c[point_y]);
								}
							}
						}
					}
				}
			}
		}
	}

	//Correlator, 2 by Def, 2 by T_12+T_21  3 by T_12+T_13+T_23 -->/12,
	//Volume for y /VOL4D, /2/pi^2*L_z^2 for derivation, *beta^2 / Nc^2 for T_munu, *2 for y_3<x_3

	hmc_float norm = (hmc_float) (NSPACE * NSPACE) / (hmc_float) (VOL4D) / (hmc_float) (NC * NC) / 12.0 / PI / PI  * beta * beta;

	*kappa_karsch_val = norm * result;
}



__kernel void kappa_clover_gpu (__global Matrixsu3StorageType * gaugefield, const hmc_float beta,  __global hmc_float * kappa_clover_val)
{

	//Energy-momentum-tensor in clover-discretization
	hmc_float t_12 [VOL4D];
	hmc_float t_13 [VOL4D];
	hmc_float t_23 [VOL4D];

	for (int t = 0; t < NTIME; t++) {
		for (int n = 0; n < VOLSPACE; n++) {
			//Compute required plaquettes
			Matrix3x3 Q_22;
			Q_22 = local_Q_plaquette(gaugefield, n, t, 2, 2);
			Matrix3x3 Q_10;
			Q_10 = local_Q_plaquette(gaugefield, n, t, 1, 0);
			Matrix3x3 Q_20;
			Q_20 = local_Q_plaquette(gaugefield, n, t, 2, 0);
			Matrix3x3 Q_02;
			Q_02 = adjoint_matrix3x3 (Q_20);
			Matrix3x3 Q_21;
			Q_21 = local_Q_plaquette(gaugefield, n, t, 2, 1);
			Matrix3x3 Q_12;
			Q_12 = adjoint_matrix3x3 (Q_21);
			Matrix3x3 Q_03;
			Q_03 = local_Q_plaquette(gaugefield, n, t, 0, 3);
			Matrix3x3 Q_30;
			Q_30 = adjoint_matrix3x3 (Q_03);
			Matrix3x3 Q_13;
			Q_13 = local_Q_plaquette( gaugefield, n, t, 1, 3);
			Matrix3x3 Q_31;
			Q_31 = adjoint_matrix3x3 (Q_13);
			Matrix3x3 Q_23;
			Q_23 = local_Q_plaquette(gaugefield, n, t, 2, 3);
			Matrix3x3 Q_32;
			Q_32 = adjoint_matrix3x3 (Q_23);
			Matrix3x3 Q_11;
			Q_11 = local_Q_plaquette(gaugefield, n, t, 1, 1);

			int point = n + VOLSPACE * t;

			Matrix3x3 tmp;
			hmc_complex tmp_cmp;

			//T_12
			//alpha=0
			tmp = subtract_matrix3x3 (Q_20, Q_02);
			tmp = multiply_matrix3x3 (Q_10, tmp);
			tmp_cmp = trace_matrix3x3(tmp);
			t_12 [point] = tmp_cmp.re;
			//alpha=1
			tmp = subtract_matrix3x3 (Q_21, Q_12);
			tmp = multiply_matrix3x3 (Q_11, tmp);
			tmp_cmp = trace_matrix3x3(tmp);
			t_12 [point] += tmp_cmp.re;
			//alpha=2, vanishes
//        subtract_3x3matrix (tmp, Q_22, Q_22);
//        multiply_3x3matrix (tmp, Q_12, tmp);
//        trace_3x3matrix (tmp_cmp, tmp);
//        t_12 [point] += tmp_cmp.re;
			//alpha=3
			tmp = subtract_matrix3x3 (Q_23, Q_32);
			tmp = multiply_matrix3x3 (Q_13, tmp);
			tmp_cmp = trace_matrix3x3(tmp);
			t_12 [point] += tmp_cmp.re;

			//T_13
			//alpha=0
			tmp = subtract_matrix3x3 (Q_30, Q_03);
			tmp = multiply_matrix3x3 (Q_10, tmp);
			tmp_cmp = trace_matrix3x3(tmp);
			t_13 [point] = tmp_cmp.re;
			//alpha=1
			tmp = subtract_matrix3x3 (Q_31, Q_13);
			tmp = multiply_matrix3x3 (Q_11, tmp);
			tmp_cmp = trace_matrix3x3(tmp);
			t_13 [point] += tmp_cmp.re;
			//alpha=2
			tmp = subtract_matrix3x3 (Q_32, Q_23);
			tmp = multiply_matrix3x3 (Q_12, tmp);
			tmp_cmp = trace_matrix3x3(tmp);
			t_13 [point] += tmp_cmp.re;
			//alpha=3, vanishes
//        subtract_3x3matrix (tmp, Q_33, Q_33);
//        multiply_3x3matrix (tmp, Q_13, tmp);
//        trace_3x3matrix (tmp_cmp, tmp);
//        t_13 [point] += tmp_cmp.re;

			//T_23
			//alpha=0
			tmp = subtract_matrix3x3 (Q_30, Q_03);
			tmp = multiply_matrix3x3 (Q_20, tmp);
			tmp_cmp = trace_matrix3x3(tmp);
			t_23 [point] = tmp_cmp.re;
			//alpha=1
			tmp = subtract_matrix3x3 (Q_31, Q_13);
			tmp = multiply_matrix3x3 (Q_21, tmp);
			tmp_cmp = trace_matrix3x3(tmp);
			t_23 [point] += tmp_cmp.re;
			//alpha=2
			tmp = subtract_matrix3x3 (Q_32, Q_23);
			tmp = multiply_matrix3x3 (Q_22, tmp);
			tmp_cmp = trace_matrix3x3(tmp);
			t_23 [point] += tmp_cmp.re;
			//alpha=3, vanishes
//        subtract_3x3matrix (tmp, Q_33, Q_33);
//        multiply_3x3matrix (tmp, Q_23, tmp);
//        trace_3x3matrix (tmp_cmp, tmp);
//        t_23 [point] += tmp_cmp.re;

		}
	}

	//Momentum
	const hmc_float deltak = 2.0 * PI / (hmc_float) NSPACE;
	hmc_float result = 0.0;


	for (int x_3 = 0; x_3 < NSPACE; x_3++) {
		for (int y_3 = 0; y_3 < x_3; y_3++) {
			hmc_float factor = 1.0 - cos(deltak * (hmc_float) (x_3 - y_3));
			for (int x_t = 0; x_t < NTIME; x_t++) {
				for (int y_t = 0; y_t < NTIME; y_t++) {
					uint3 coord_y, coord_y;
					for (int coord_x.x = 0; coord_x.x < NSPACE; coord_x.x++) {
						for (int coord_y.x = 0; coord_y.x < NSPACE; coord_y.x++) {
							for (int coord_x.y = 0; coord_x.y < NSPACE; coord_x.y++) {
								for (int coord_y.y = 0; coord_y.y < NSPACE; coord_y.y++) {

									//new method get_tnspace which gives n_x+VOLSPACE*x_t
									int n_x = get_nspace (coord_x);
									int point_x = n_x + VOLSPACE * x_t;
									int n_y = get_nspace (coord_y);
									int point_y = n_y + VOLSPACE * y_t;

									//(T_12(x) T_12(y) + T_21(x) T_21(y) + T_13(x) T_13(y)) * factor
									result += factor  * ( t_12[point_x] * t_12[point_y]
									                      + t_13[point_x] * t_13[point_y]
									                      + t_23[point_x] * t_23[point_y]);
								}
							}
						}
					}
				}
			}
		}
	}

	//Normalization
	// 1/3 for averaging T_ij, 1/V/Nt for averaging y, L^2/2/pi^2 for derivation, (-1/64)^2 for Clover and T_munu^2, beta^2/Nc^2 for T_munu^2
	// *2 for temp + conj (temp) *2 for for-loop
	// = beta^2 * L^2/ (55296 * V * Nt * pi^2)
	hmc_float norm = (hmc_float) (NSPACE * NSPACE) / (hmc_float) (VOL4D) / PI / PI * beta * beta / 55296. ;

	* kappa_clover_val = norm * result;


}
