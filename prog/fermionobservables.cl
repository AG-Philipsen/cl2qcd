/**
 @file fermion-observables
*/

// //this is the pseudoscalar pion correlator in z-direction from pointsources
__kernel void correlator_ps_z(__global const spinor * const restrict phi, __global hmc_float * const restrict out)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NSPACE; id_tmp += global_size) {
		hmc_float correlator = 0.;
		uint3 coord;
		coord.z = id_tmp;
		for(int k = 0; k < NUM_SOURCES; k++) {
			for(int t = 0; t < NTIME; t++) {
				for(coord.x = 0; coord.x < NSPACE; coord.x++) {
					for(coord.y = 0; coord.y < NSPACE; coord.y++) {
						int nspace = get_nspace(coord);
						spinor tmp = phi[get_global_pos(nspace, t) + VOL4D * k];
						correlator += spinor_squarenorm(tmp);
					}
				}
			}
		}
		//now, this should finally be the correct normalisation for the physical fields
		//one factor of 2*kappa per field and we construct the correlator from a multiplication of two fields phi

		hmc_float fac = NSPACE * NSPACE * NTIME;
		out[id_tmp] = 2. * KAPPA * 2.* KAPPA * correlator / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int z=0; z<NSPACE; z++)
	//  printf("%i\t(%.12e)\n", z, out[z]);
	//     }
	//#endif

}

// //this is the pseudoscalar pion correlator in t-direction from pointsources
__kernel void correlator_ps_t(__global const spinor * const restrict phi, __global hmc_float * const restrict out)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NTIME threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NTIME; id_tmp += global_size) {
		hmc_float correlator = 0.;
		uint3 coord;
		int t = id_tmp;
		for(int k = 0; k < NUM_SOURCES; k++) {
			for(coord.z = 0; coord.z < NSPACE; coord.z++) {
				for(coord.x = 0; coord.x < NSPACE; coord.x++) {
					for(coord.y = 0; coord.y < NSPACE; coord.y++) {
						int nspace = get_nspace(coord);
						spinor tmp = phi[get_global_pos(nspace, t) + VOL4D * k];
						correlator += spinor_squarenorm(tmp);
					}
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NSPACE;
		out[id_tmp] = 2. * KAPPA * 2. * KAPPA * correlator / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int t=0; t<NTIME; t++)
	//  printf("%i\t(%.12e)\n", t, out[t]);
	//     }
	//#endif

}


// //this is the scalar correlator in z-direction from pointsources
__kernel void correlator_sc_z(__global const spinor * const restrict phi, __global hmc_float * const restrict out)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NSPACE; id_tmp += global_size) {
		hmc_float correlator = 0.;
		uint3 coord;
		coord.z = id_tmp;
		for(int t = 0; t < NTIME; t++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					int alpha;
					int k;
					spinor tmp;

					for(int color = 0; color < 3; color++) {

						alpha = 0;
						k = spinor_element(alpha, color);
						tmp = phi[get_global_pos(nspace, t) + VOL4D * k];
						correlator += - su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);

						alpha = 1;
						k = spinor_element(alpha, color);
						tmp = phi[get_global_pos(nspace, t) + VOL4D * k];
						correlator += - su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);

						alpha = 2;
						k = spinor_element(alpha, color);
						tmp = phi[get_global_pos(nspace, t) + VOL4D * k];
						correlator += su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);

						alpha = 3;
						k = spinor_element(alpha, color);
						tmp = phi[get_global_pos(nspace, t) + VOL4D * k];
						correlator += su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);
					}
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NTIME;
		out[id_tmp] = 2. * KAPPA * 2. * KAPPA * correlator / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int z=0; z<NSPACE; z++)
	//  printf("%i\t(%.12e)\n", z, out[z]);
	//     }
	//#endif

}

// //this is the scalar correlator in t-direction from pointsources
__kernel void correlator_sc_t(__global const spinor * const restrict phi, __global hmc_float * const restrict out)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NTIME threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NTIME; id_tmp += global_size) {
		hmc_float correlator = 0.;
		uint3 coord;
		int t = id_tmp;
		for(coord.z = 0; coord.z < NSPACE; coord.z++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					int alpha;
					int k;
					spinor tmp;

					for(int color = 0; color < 3; color++) {

						alpha = 0;
						k = spinor_element(alpha, color);
						tmp = phi[get_global_pos(nspace, t) + VOL4D * k];
						correlator += - su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);

						alpha = 1;
						k = spinor_element(alpha, color);
						tmp = phi[get_global_pos(nspace, t) + VOL4D * k];
						correlator += - su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);

						alpha = 2;
						k = spinor_element(alpha, color);
						tmp = phi[get_global_pos(nspace, t) + VOL4D * k];
						correlator += su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);

						alpha = 3;
						k = spinor_element(alpha, color);
						tmp = phi[get_global_pos(nspace, t) + VOL4D * k];
						correlator += su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);
					}
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NSPACE;
		out[id_tmp] = 2. * KAPPA * 2.*KAPPA * correlator / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int t=0; t<NTIME; t++)
	//  printf("%i\t(%.12e)\n", t, out[t]);
	//     }
	//endif

}


// //this is the vector correlator in z-direction from pointsources (gamma4*gamma1)
__kernel void correlator_vx_z(__global const spinor * const restrict phi, __global hmc_float * const restrict out)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NSPACE; id_tmp += global_size) {
		hmc_complex correlator;
		correlator.re = 0.0f;
		correlator.im = 0.0f;
		uint3 coord;
		coord.z = id_tmp;

		for(int t = 0; t < NTIME; t++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					int alpha;
					int k;
					spinor tmp_a;
					spinor tmp_b;
					hmc_complex restmp;

					for(int color = 0; color < 3; color++) {

						alpha = 0;
						k = spinor_element(alpha, color);
						tmp_a = phi[get_global_pos(nspace, t) + VOL4D * k];
						alpha = 1;
						k = spinor_element(alpha, color);
						tmp_b = phi[get_global_pos(nspace, t) + VOL4D * k];

						restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
						correlator.re += restmp.re;
						correlator.im += restmp.im;


						alpha = 2;
						k = spinor_element(alpha, color);
						tmp_a = phi[get_global_pos(nspace, t) + VOL4D * k];
						alpha = 3;
						k = spinor_element(alpha, color);
						tmp_b = phi[get_global_pos(nspace, t) + VOL4D * k];

						restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
						correlator.re += restmp.re;
						correlator.im += restmp.im;
					}
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NTIME;
		out[id_tmp] = 2. * KAPPA * 2. * KAPPA * 2.*correlator.re / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int z=0; z<NSPACE; z++)
	//  printf("%i\t(%.12e)\n", z, out[z]);
	//     }
	//#endif

}

// //this is the vector correlator in t-direction from pointsources (gamma4*gamma1)
__kernel void correlator_vx_t(__global const spinor * const restrict phi, __global hmc_float * const restrict out)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NTIME; id_tmp += global_size) {
		hmc_complex correlator;
		correlator.re = 0.0f;
		correlator.im = 0.0f;
		uint3 coord;
		int t = id_tmp;
		for(coord.z = 0; coord.z < NSPACE; coord.z++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					int alpha;
					int k;
					spinor tmp_a;
					spinor tmp_b;
					hmc_complex restmp;

					for(int color = 0; color < 3; color++) {

						alpha = 0;
						k = spinor_element(alpha, color);
						tmp_a = phi[get_global_pos(nspace, t) + VOL4D * k];
						alpha = 1;
						k = spinor_element(alpha, color);
						tmp_b = phi[get_global_pos(nspace, t) + VOL4D * k];

						restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
						correlator.re += restmp.re;
						correlator.im += restmp.im;


						alpha = 2;
						k = spinor_element(alpha, color);
						tmp_a = phi[get_global_pos(nspace, t) + VOL4D * k];
						alpha = 3;
						k = spinor_element(alpha, color);
						tmp_b = phi[get_global_pos(nspace, t) + VOL4D * k];

						restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
						correlator.re += restmp.re;
						correlator.im += restmp.im;
					}
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NSPACE;
		out[id_tmp] = 2. * KAPPA * 2. * KAPPA * 2. * correlator.re / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int t=0; t<NTIME; t++)
	//  printf("%i\t(%.12e)\n", t, out[t]);
	//     }
	//#endif

}


// //this is the vector correlator in z-direction from pointsources (gamma4*gamma2)
__kernel void correlator_vy_z(__global const spinor * const restrict phi, __global hmc_float * const restrict out)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NSPACE; id_tmp += global_size) {
		hmc_complex correlator;
		correlator.re = 0.0f;
		correlator.im = 0.0f;
		uint3 coord;
		coord.z = id_tmp;
		for(int t = 0; t < NTIME; t++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					int alpha;
					int k;
					spinor tmp_a;
					spinor tmp_b;
					hmc_complex restmp;

					for(int color = 0; color < 3; color++) {

						alpha = 0;
						k = spinor_element(alpha, color);
						tmp_a = phi[get_global_pos(nspace, t) + VOL4D * k];
						alpha = 1;
						k = spinor_element(alpha, color);
						tmp_b = phi[get_global_pos(nspace, t) + VOL4D * k];

						restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
						correlator.re -= restmp.re;
						correlator.im -= restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
						correlator.re -= restmp.re;
						correlator.im -= restmp.im;


						alpha = 2;
						k = spinor_element(alpha, color);
						tmp_a = phi[get_global_pos(nspace, t) + VOL4D * k];
						alpha = 3;
						k = spinor_element(alpha, color);
						tmp_b = phi[get_global_pos(nspace, t) + VOL4D * k];

						restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
						correlator.re -= restmp.re;
						correlator.im -= restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
						correlator.re -= restmp.re;
						correlator.im -= restmp.im;
					}
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NTIME;
		out[id_tmp] = 2. * KAPPA * 2. * KAPPA * 2.*correlator.re / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int z=0; z<NSPACE; z++)
	//  printf("%i\t(%.12e)\n", z, out[z]);
	//     }
	//#endif

}

// //this is the vector correlator in t-direction from pointsources (gamma4*gamma2)
__kernel void correlator_vy_t(__global const spinor * const restrict phi, __global hmc_float * const restrict out)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NTIME; id_tmp += global_size) {
		hmc_complex correlator;
		correlator.re = 0.0f;
		correlator.im = 0.0f;
		uint3 coord;
		int t = id_tmp;
		for(coord.z = 0; coord.z < NSPACE; coord.z++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					int alpha;
					int k;
					spinor tmp_a;
					spinor tmp_b;
					hmc_complex restmp;

					for(int color = 0; color < 3; color++) {

						alpha = 0;
						k = spinor_element(alpha, color);
						tmp_a = phi[get_global_pos(nspace, t) + VOL4D * k];
						alpha = 1;
						k = spinor_element(alpha, color);
						tmp_b = phi[get_global_pos(nspace, t) + VOL4D * k];

						restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
						correlator.re -= restmp.re;
						correlator.im -= restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
						correlator.re -= restmp.re;
						correlator.im -= restmp.im;


						alpha = 2;
						k = spinor_element(alpha, color);
						tmp_a = phi[get_global_pos(nspace, t) + VOL4D * k];
						alpha = 3;
						k = spinor_element(alpha, color);
						tmp_b = phi[get_global_pos(nspace, t) + VOL4D * k];

						restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
						correlator.re -= restmp.re;
						correlator.im -= restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
						correlator.re -= restmp.re;
						correlator.im -= restmp.im;
					}
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NSPACE;
		out[id_tmp] = 2. * KAPPA * 2. * KAPPA * 2. * correlator.re / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int t=0; t<NTIME; t++)
	//  printf("%i\t(%.12e)\n", t, out[t]);
	//     }
	//#endif

}



// //this is the vector correlator in z-direction from pointsources (gamma4*gamma3)
__kernel void correlator_vz_z(__global const spinor * const restrict phi, __global hmc_float * const restrict out)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NSPACE; id_tmp += global_size) {
		hmc_float correlator = 0.0f;
		uint3 coord;
		coord.z = id_tmp;
		for(int t = 0; t < NTIME; t++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					int alpha;
					int k;
					spinor tmp;

					for(int color = 0; color < 3; color++) {

						alpha = 0;
						k = spinor_element(alpha, color);
						tmp = phi[get_global_pos(nspace, t) + VOL4D * k];
						correlator += su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);

						alpha = 1;
						k = spinor_element(alpha, color);
						tmp = phi[get_global_pos(nspace, t) + VOL4D * k];
						correlator += -su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);

						alpha = 2;
						k = spinor_element(alpha, color);
						tmp = phi[get_global_pos(nspace, t) + VOL4D * k];
						correlator += su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);

						alpha = 3;
						k = spinor_element(alpha, color);
						tmp = phi[get_global_pos(nspace, t) + VOL4D * k];
						correlator += -su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);
					}
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NTIME;
		out[id_tmp] = 2. * KAPPA * 2.*KAPPA * correlator / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int z=0; z<NSPACE; z++)
	//  printf("%i\t(%.12e)\n", z, out[z]);
	//     }
	//#endif

}

// //this is the vector correlator in t-direction from pointsources (gamma4*gamma3)
__kernel void correlator_vz_t(__global const spinor * const restrict phi, __global hmc_float * const restrict out)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NTIME; id_tmp += global_size) {
		hmc_float correlator = 0.0f;
		uint3 coord;
		int t = id_tmp;
		for(coord.z = 0; coord.z < NSPACE; coord.z++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					int alpha;
					int k;
					spinor tmp;

					for(int color = 0; color < 3; color++) {

						alpha = 0;
						k = spinor_element(alpha, color);
						tmp = phi[get_global_pos(nspace, t) + VOL4D * k];
						correlator += su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);

						alpha = 1;
						k = spinor_element(alpha, color);
						tmp = phi[get_global_pos(nspace, t) + VOL4D * k];
						correlator += -su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);

						alpha = 2;
						k = spinor_element(alpha, color);
						tmp = phi[get_global_pos(nspace, t) + VOL4D * k];
						correlator += su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);

						alpha = 3;
						k = spinor_element(alpha, color);
						tmp = phi[get_global_pos(nspace, t) + VOL4D * k];
						correlator += -su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);
					}
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NSPACE;
		out[id_tmp] = 2. * KAPPA * 2. * KAPPA * correlator / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int t=0; t<NTIME; t++)
	//  printf("%i\t(%.12e)\n", t, out[t]);
	//     }
	//#endif

}

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// axial vector correlators

// //this is the axial vector correlator in z-direction from pointsources (gamma5*gamma4*gamma1)
__kernel void correlator_ax_z(__global const spinor * const restrict phi, __global hmc_float * const restrict out)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NSPACE; id_tmp += global_size) {
		hmc_complex correlator;
		correlator.re = 0.0f;
		correlator.im = 0.0f;
		uint3 coord;
		coord.z = id_tmp;
		for(int t = 0; t < NTIME; t++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					int alpha;
					int k;
					spinor tmp_a;
					spinor tmp_b;
					hmc_complex restmp;

					for(int color = 0; color < 3; color++) {

						alpha = 0;
						k = spinor_element(alpha, color);
						tmp_a = phi[get_global_pos(nspace, t) + VOL4D * k];
						alpha = 1;
						k = spinor_element(alpha, color);
						tmp_b = phi[get_global_pos(nspace, t) + VOL4D * k];

						restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
						correlator.re -= restmp.re;
						correlator.im -= restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
						correlator.re -= restmp.re;
						correlator.im -= restmp.im;


						alpha = 2;
						k = spinor_element(alpha, color);
						tmp_a = phi[get_global_pos(nspace, t) + VOL4D * k];
						alpha = 3;
						k = spinor_element(alpha, color);
						tmp_b = phi[get_global_pos(nspace, t) + VOL4D * k];

						restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
						correlator.re -= restmp.re;
						correlator.im -= restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
						correlator.re -= restmp.re;
						correlator.im -= restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
						correlator.re += restmp.re;
						correlator.im += restmp.im;
					}
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NTIME;
		out[id_tmp] = - 2. * KAPPA * 2.*KAPPA * 2.*correlator.re / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int z=0; z<NSPACE; z++)
	//  printf("%i\t(%.12e)\n", z, out[z]);
	//     }
	//#endif

}

// //this is the axial vector correlator in t-direction from pointsources (gamma5*gamma4*gamma1)
__kernel void correlator_ax_t(__global const spinor * const restrict phi, __global hmc_float * const restrict out)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NTIME; id_tmp += global_size) {
		hmc_complex correlator;
		correlator.re = 0.0f;
		correlator.im = 0.0f;
		uint3 coord;
		int t = id_tmp;
		for(coord.z = 0; coord.z < NSPACE; coord.z++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					int alpha;
					int k;
					spinor tmp_a;
					spinor tmp_b;
					hmc_complex restmp;

					for(int color = 0; color < 3; color++) {

						alpha = 0;
						k = spinor_element(alpha, color);
						tmp_a = phi[get_global_pos(nspace, t) + VOL4D * k];
						alpha = 1;
						k = spinor_element(alpha, color);
						tmp_b = phi[get_global_pos(nspace, t) + VOL4D * k];

						restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
						correlator.re -= restmp.re;
						correlator.im -= restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
						correlator.re -= restmp.re;
						correlator.im -= restmp.im;


						alpha = 2;
						k = spinor_element(alpha, color);
						tmp_a = phi[get_global_pos(nspace, t) + VOL4D * k];
						alpha = 3;
						k = spinor_element(alpha, color);
						tmp_b = phi[get_global_pos(nspace, t) + VOL4D * k];

						restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
						correlator.re -= restmp.re;
						correlator.im -= restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
						correlator.re -= restmp.re;
						correlator.im -= restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
						correlator.re += restmp.re;
						correlator.im += restmp.im;
					}
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NSPACE;
		out[id_tmp] = - 2. * KAPPA * 2.*KAPPA * 2. * correlator.re / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int t=0; t<NTIME; t++)
	//  printf("%i\t(%.12e)\n", t, out[t]);
	//     }
	//#endif

}


// //this is the axial vector correlator in z-direction from pointsources (gamma5*gamma4*gamma2)
__kernel void correlator_ay_z(__global const spinor * const restrict phi, __global hmc_float * const restrict out)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NSPACE; id_tmp += global_size) {
		hmc_complex correlator;
		correlator.re = 0.0f;
		correlator.im = 0.0f;
		uint3 coord;
		coord.z = id_tmp;
		for(int t = 0; t < NTIME; t++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					int alpha;
					int k;
					spinor tmp_a;
					spinor tmp_b;
					hmc_complex restmp;

					for(int color = 0; color < 3; color++) {

						alpha = 0;
						k = spinor_element(alpha, color);
						tmp_a = phi[get_global_pos(nspace, t) + VOL4D * k];
						alpha = 1;
						k = spinor_element(alpha, color);
						tmp_b = phi[get_global_pos(nspace, t) + VOL4D * k];

						restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
						correlator.re -= restmp.re;
						correlator.im -= restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
						correlator.re -= restmp.re;
						correlator.im -= restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
						correlator.re += restmp.re;
						correlator.im += restmp.im;


						alpha = 2;
						k = spinor_element(alpha, color);
						tmp_a = phi[get_global_pos(nspace, t) + VOL4D * k];
						alpha = 3;
						k = spinor_element(alpha, color);
						tmp_b = phi[get_global_pos(nspace, t) + VOL4D * k];

						restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
						correlator.re -= restmp.re;
						correlator.im -= restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
						correlator.re -= restmp.re;
						correlator.im -= restmp.im;
					}
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NTIME;
		out[id_tmp] = - 2. * KAPPA * 2.*KAPPA * 2.*correlator.re / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int z=0; z<NSPACE; z++)
	//  printf("%i\t(%.12e)\n", z, out[z]);
	//     }
	//#endif

}

// //this is the axial vector correlator in t-direction from pointsources (gamma5*gamma4*gamma2)
__kernel void correlator_ay_t(__global const spinor * const restrict phi, __global hmc_float * const restrict out)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NTIME; id_tmp += global_size) {
		hmc_complex correlator;
		correlator.re = 0.0f;
		correlator.im = 0.0f;
		uint3 coord;
		int t = id_tmp;
		for(coord.z = 0; coord.z < NSPACE; coord.z++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					int alpha;
					int k;
					spinor tmp_a;
					spinor tmp_b;
					hmc_complex restmp;

					for(int color = 0; color < 3; color++) {

						alpha = 0;
						k = spinor_element(alpha, color);
						tmp_a = phi[get_global_pos(nspace, t) + VOL4D * k];
						alpha = 1;
						k = spinor_element(alpha, color);
						tmp_b = phi[get_global_pos(nspace, t) + VOL4D * k];

						restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
						correlator.re -= restmp.re;
						correlator.im -= restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
						correlator.re -= restmp.re;
						correlator.im -= restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
						correlator.re += restmp.re;
						correlator.im += restmp.im;


						alpha = 2;
						k = spinor_element(alpha, color);
						tmp_a = phi[get_global_pos(nspace, t) + VOL4D * k];
						alpha = 3;
						k = spinor_element(alpha, color);
						tmp_b = phi[get_global_pos(nspace, t) + VOL4D * k];

						restmp = su3vec_scalarproduct(tmp_a.e0, tmp_b.e1);
						correlator.re -= restmp.re;
						correlator.im -= restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e1, tmp_b.e0);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e2, tmp_b.e3);
						correlator.re += restmp.re;
						correlator.im += restmp.im;

						restmp = su3vec_scalarproduct(tmp_a.e3, tmp_b.e2);
						correlator.re -= restmp.re;
						correlator.im -= restmp.im;
					}
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NSPACE;
		out[id_tmp] = - 2. * KAPPA * 2.*KAPPA * 2. * correlator.re / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int t=0; t<NTIME; t++)
	//  printf("%i\t(%.12e)\n", t, out[t]);
	//     }
	//#endif

}



// //this is the axial vector correlator in z-direction from pointsources (gamma5*gamma4*gamma3)
__kernel void correlator_az_z(__global const spinor * const restrict phi, __global hmc_float * const restrict out)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NSPACE; id_tmp += global_size) {
		hmc_float correlator = 0.0f;
		uint3 coord;
		coord.z = id_tmp;
		for(int t = 0; t < NTIME; t++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					int alpha;
					int k;
					spinor tmp;

					for(int color = 0; color < 3; color++) {

						alpha = 0;
						k = spinor_element(alpha, color);
						tmp = phi[get_global_pos(nspace, t) + VOL4D * k];
						correlator += su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);

						alpha = 1;
						k = spinor_element(alpha, color);
						tmp = phi[get_global_pos(nspace, t) + VOL4D * k];
						correlator += -su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);

						alpha = 2;
						k = spinor_element(alpha, color);
						tmp = phi[get_global_pos(nspace, t) + VOL4D * k];
						correlator += -su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);

						alpha = 3;
						k = spinor_element(alpha, color);
						tmp = phi[get_global_pos(nspace, t) + VOL4D * k];
						correlator += su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);
					}
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NTIME;
		out[id_tmp] = - 2. * KAPPA * 2.*KAPPA * correlator / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int z=0; z<NSPACE; z++)
	//  printf("%i\t(%.12e)\n", z, out[z]);
	//     }
	//#endif

}

// //this is the axial vector correlator in t-direction from pointsources (gamma5*gamma4*gamma3)
__kernel void correlator_az_t(__global const spinor * const restrict phi, __global hmc_float * const restrict out)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NTIME threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NTIME; id_tmp += global_size) {
		hmc_float correlator = 0.0f;
		uint3 coord;
		int t = id_tmp;
		for(coord.z = 0; coord.z < NSPACE; coord.z++) {
			for(coord.x = 0; coord.x < NSPACE; coord.x++) {
				for(coord.y = 0; coord.y < NSPACE; coord.y++) {
					int nspace = get_nspace(coord);
					int alpha;
					int k;
					spinor tmp;

					for(int color = 0; color < 3; color++) {

						alpha = 0;
						k = spinor_element(alpha, color);
						tmp = phi[get_global_pos(nspace, t) + VOL4D * k];
						correlator += su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);

						alpha = 1;
						k = spinor_element(alpha, color);
						tmp = phi[get_global_pos(nspace, t) + VOL4D * k];
						correlator += -su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);

						alpha = 2;
						k = spinor_element(alpha, color);
						tmp = phi[get_global_pos(nspace, t) + VOL4D * k];
						correlator += -su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);

						alpha = 3;
						k = spinor_element(alpha, color);
						tmp = phi[get_global_pos(nspace, t) + VOL4D * k];
						correlator += su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);
					}
				}
			}
		}
		hmc_float fac = NSPACE * NSPACE * NSPACE;
		out[id_tmp] = - 2. * KAPPA * 2.* KAPPA * correlator / fac;
	}


	//LZ: print directly to stdout for debugging:
	//#ifdef ENABLE_PRINTF
	//     if(id == 0) {
	//       for(int t=0; t<NTIME; t++)
	//  printf("%i\t(%.12e)\n", t, out[t]);
	//     }
	//#endif

}
