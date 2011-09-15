/**
 @file fermion-observables
*/

// //this is the pseudoscalar pion correlator in z-direction from pointsources
__kernel void correlator_ps_z(__global spinorfield* phi, __global hmc_float * out)
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
		int coord[4]; //LZ: int4 would be nicer but that cannot go into the current get_nspace() function...
		coord[3] = id_tmp;
		for(int k = 0; k < NUM_SOURCES; k++) {
			for(int t = 0; t < NTIME; t++) {
				for(int x = 0; x < NSPACE; x++) {
					for(int y = 0; y < NSPACE; y++) {
						coord[1] = x;
						coord[2] = y;
						int nspace = get_nspace(coord);
						spinor tmp = phi[get_global_pos(nspace, t) + VOL4D*k];
						correlator += spinor_squarenorm(tmp);
					}
				}
			}
		}
		out[id_tmp] = 4.*KAPPA * KAPPA * correlator;
	}


	//LZ: print directly to stdout for debugging:
	//     if(id == 0) {
	//       for(int z=0; z<NSPACE; z++)
	//  printf("%i\t(%.12e)\n", z, out[z]);
	//     }

}

// //this is the pseudoscalar pion correlator in t-direction from pointsources
__kernel void correlator_ps_t(__global spinorfield* phi, __global hmc_float * out)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NTIME; id_tmp += global_size) {
		hmc_float correlator = 0.;
		int coord[4]; //LZ: int4 would be nicer but that cannot go into the current get_nspace() function...
		int t = id_tmp;
		for(int k = 0; k < NUM_SOURCES; k++) {
			for(int z = 0; z < NSPACE; z++) {
				for(int x = 0; x < NSPACE; x++) {
					for(int y = 0; y < NSPACE; y++) {
						coord[1] = x;
						coord[2] = y;
						coord[3] = z;
						int nspace = get_nspace(coord);
						spinor tmp = phi[get_global_pos(nspace, t) + VOL4D*k];
						correlator += spinor_squarenorm(tmp);
					}
				}
			}
		}
		out[id_tmp] = 4.*KAPPA * KAPPA * correlator;
	}


	//LZ: print directly to stdout for debugging:
	//     if(id == 0) {
	//       for(int t=0; t<NTIME; t++)
	//  printf("%i\t(%.12e)\n", t, out[t]);
	//     }

}


// //this is the scalar correlator in z-direction from pointsources
__kernel void correlator_sc_z(__global spinorfield* phi, __global hmc_float * out)
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
		int coord[4]; //LZ: int4 would be nicer but that cannot go into the current get_nspace() function...
		coord[3] = id_tmp;
		for(int t = 0; t < NTIME; t++) {
			for(int x = 0; x < NSPACE; x++) {
				for(int y = 0; y < NSPACE; y++) {
					coord[1] = x;
					coord[2] = y;
					int nspace = get_nspace(coord);
					int k;
					spinor tmp;

					k = 0;
					tmp = phi[get_global_pos(nspace, t) + VOL4D*k];
					correlator += - su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);

					k = 1;
					tmp = phi[get_global_pos(nspace, t) + VOL4D*k];
					correlator += - su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);

					k = 2;
					tmp = phi[get_global_pos(nspace, t) + VOL4D*k];
					correlator += su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);

					k = 2;
					tmp = phi[get_global_pos(nspace, t) + VOL4D*k];
					correlator += su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);

				}
			}
		}
		out[id_tmp] = 4.*KAPPA * KAPPA * correlator;
	}


	//LZ: print directly to stdout for debugging:
	//     if(id == 0) {
	//       for(int z=0; z<NSPACE; z++)
	//  printf("%i\t(%.12e)\n", z, out[z]);
	//     }

}

// //this is the scalar correlator in t-direction from pointsources
__kernel void correlator_sc_t(__global spinorfield* phi, __global hmc_float * out)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	//suppose that there are NSPACE threads (one for each entry of the correlator)
	for(int id_tmp = id; id_tmp < NTIME; id_tmp += global_size) {
		hmc_float correlator = 0.;
		int coord[4]; //LZ: int4 would be nicer but that cannot go into the current get_nspace() function...
		int t = id_tmp;
		for(int z = 0; z < NSPACE; z++) {
			for(int x = 0; x < NSPACE; x++) {
				for(int y = 0; y < NSPACE; y++) {
					coord[1] = x;
					coord[2] = y;
					coord[3] = z;
					int nspace = get_nspace(coord);
					int k;
					spinor tmp;

					k = 0;
					tmp = phi[get_global_pos(nspace, t) + VOL4D*k];
					correlator += - su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);

					k = 1;
					tmp = phi[get_global_pos(nspace, t) + VOL4D*k];
					correlator += - su3vec_squarenorm(tmp.e0) - su3vec_squarenorm(tmp.e1) + su3vec_squarenorm(tmp.e2) + su3vec_squarenorm(tmp.e3);

					k = 2;
					tmp = phi[get_global_pos(nspace, t) + VOL4D*k];
					correlator += su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);

					k = 2;
					tmp = phi[get_global_pos(nspace, t) + VOL4D*k];
					correlator += su3vec_squarenorm(tmp.e0) + su3vec_squarenorm(tmp.e1) - su3vec_squarenorm(tmp.e2) - su3vec_squarenorm(tmp.e3);

				}
			}
		}
		out[id_tmp] = 4.*KAPPA * KAPPA * correlator;
	}


	//LZ: print directly to stdout for debugging:
	//     if(id == 0) {
	//       for(int t=0; t<NTIME; t++)
	//  printf("%i\t(%.12e)\n", t, out[t]);
	//     }

}


// //this is the vector correlator in z-direction from pointsources (gamma4*gamma1)
__kernel void correlator_vx_z(__global spinorfield* phi, __global hmc_float * out)
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
		int coord[4]; //LZ: int4 would be nicer but that cannot go into the current get_nspace() function...
		coord[3] = id_tmp;
		for(int t = 0; t < NTIME; t++) {
			for(int x = 0; x < NSPACE; x++) {
				for(int y = 0; y < NSPACE; y++) {
					coord[1] = x;
					coord[2] = y;
					int nspace = get_nspace(coord);
					int k;
					spinor tmp_a;
					spinor tmp_b;
					hmc_complex restmp;

					k = 0;
					tmp_a = phi[get_global_pos(nspace, t) + VOL4D*k];
					k = 1;
					tmp_b = phi[get_global_pos(nspace, t) + VOL4D*k];

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


					k = 2;
					tmp_a = phi[get_global_pos(nspace, t) + VOL4D*k];
					k = 3;
					tmp_b = phi[get_global_pos(nspace, t) + VOL4D*k];

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
		out[id_tmp] = 4.*KAPPA * KAPPA * 2.*correlator.re;
	}


	//LZ: print directly to stdout for debugging:
	//     if(id == 0) {
	//       for(int z=0; z<NSPACE; z++)
	//  printf("%i\t(%.12e)\n", z, out[z]);
	//     }

}

// //this is the vector correlator in t-direction from pointsources (gamma4*gamma1)
__kernel void correlator_vx_t(__global spinorfield* phi, __global hmc_float * out)
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
		int coord[4]; //LZ: int4 would be nicer but that cannot go into the current get_nspace() function...
		int t = id_tmp;
		for(int z = 0; z < NSPACE; z++) {
			for(int x = 0; x < NSPACE; x++) {
				for(int y = 0; y < NSPACE; y++) {
					coord[1] = x;
					coord[2] = y;
					coord[3] = z;
					int nspace = get_nspace(coord);
					int k;
					spinor tmp_a;
					spinor tmp_b;
					hmc_complex restmp;


					k = 0;
					tmp_a = phi[get_global_pos(nspace, t) + VOL4D*k];
					k = 1;
					tmp_b = phi[get_global_pos(nspace, t) + VOL4D*k];

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


					k = 2;
					tmp_a = phi[get_global_pos(nspace, t) + VOL4D*k];
					k = 3;
					tmp_b = phi[get_global_pos(nspace, t) + VOL4D*k];

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
		out[id_tmp] = 4.*KAPPA * KAPPA * 2. * correlator.re;
	}


	//LZ: print directly to stdout for debugging:
	//     if(id == 0) {
	//       for(int t=0; t<NTIME; t++)
	//  printf("%i\t(%.12e)\n", t, out[t]);
	//     }

}


// //this is the vector correlator in z-direction from pointsources (gamma4*gamma2)
__kernel void correlator_vy_z(__global spinorfield* phi, __global hmc_float * out)
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
		int coord[4]; //LZ: int4 would be nicer but that cannot go into the current get_nspace() function...
		coord[3] = id_tmp;
		for(int t = 0; t < NTIME; t++) {
			for(int x = 0; x < NSPACE; x++) {
				for(int y = 0; y < NSPACE; y++) {
					coord[1] = x;
					coord[2] = y;
					int nspace = get_nspace(coord);
					int k;
					spinor tmp_a;
					spinor tmp_b;
					hmc_complex restmp;

					k = 0;
					tmp_a = phi[get_global_pos(nspace, t) + VOL4D*k];
					k = 1;
					tmp_b = phi[get_global_pos(nspace, t) + VOL4D*k];

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


					k = 2;
					tmp_a = phi[get_global_pos(nspace, t) + VOL4D*k];
					k = 3;
					tmp_b = phi[get_global_pos(nspace, t) + VOL4D*k];

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
		out[id_tmp] = 4.*KAPPA * KAPPA * 2.*correlator.re;
	}


	//LZ: print directly to stdout for debugging:
	//     if(id == 0) {
	//       for(int z=0; z<NSPACE; z++)
	//  printf("%i\t(%.12e)\n", z, out[z]);
	//     }

}

// //this is the vector correlator in t-direction from pointsources (gamma4*gamma2)
__kernel void correlator_vy_t(__global spinorfield* phi, __global hmc_float * out)
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
		int coord[4]; //LZ: int4 would be nicer but that cannot go into the current get_nspace() function...
		int t = id_tmp;
		for(int z = 0; z < NSPACE; z++) {
			for(int x = 0; x < NSPACE; x++) {
				for(int y = 0; y < NSPACE; y++) {
					coord[1] = x;
					coord[2] = y;
					coord[3] = z;
					int nspace = get_nspace(coord);
					int k;
					spinor tmp_a;
					spinor tmp_b;
					hmc_complex restmp;


					k = 0;
					tmp_a = phi[get_global_pos(nspace, t) + VOL4D*k];
					k = 1;
					tmp_b = phi[get_global_pos(nspace, t) + VOL4D*k];

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


					k = 2;
					tmp_a = phi[get_global_pos(nspace, t) + VOL4D*k];
					k = 3;
					tmp_b = phi[get_global_pos(nspace, t) + VOL4D*k];

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
		out[id_tmp] = 4.*KAPPA * KAPPA * 2. * correlator.re;
	}


	//LZ: print directly to stdout for debugging:
	//     if(id == 0) {
	//       for(int t=0; t<NTIME; t++)
	//  printf("%i\t(%.12e)\n", t, out[t]);
	//     }

}

