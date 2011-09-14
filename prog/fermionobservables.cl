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
		hmc_float correlator_ps = 0.;
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
						correlator_ps += spinor_squarenorm(tmp);
					}
				}
			}
		}
		out[id_tmp] = 4.*KAPPA * KAPPA * correlator_ps;
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
		hmc_float correlator_ps = 0.;
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
						correlator_ps += spinor_squarenorm(tmp);
					}
				}
			}
		}
		out[id_tmp] = 4.*KAPPA * KAPPA * correlator_ps;
	}


	//LZ: print directly to stdout for debugging:
	//     if(id == 0) {
	//       for(int t=0; t<NTIME; t++)
	//  printf("%i\t(%.12e)\n", t, out[t]);
	//     }

}

