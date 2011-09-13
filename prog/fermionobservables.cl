/**
 @file fermion-observables
*/

//this is the pseudoscalar pion correlator in z-direction
//__kernel void ps_correlator(__global spinorfield* phi){
// int local_size = get_local_size(0);
//   int global_size = get_global_size(0);
//   int id = get_global_id(0);
//   int loc_idx = get_local_id(0);
//   int num_groups = get_num_groups(0);
//   int group_id = get_group_id (0);
//
//   if(id==0){
//     hmc_float correlator_ps[NSPACE];
//     for(int i = 0; i<NSPACE; i++){
//       correlator_ps[i]=0.;
//     }
//     for(int k=0; k<12; k++) {
//       for(int timepos = 0; timepos<NTIME; timepos++) {
//	 for(int spacepos = 0; spacepos<VOLSPACE; spacepos++) {
//		//correlator_ps[z] += |phi(n,t)|^2
//     spinor tmp = phi[get_global_pos(spacepos, timepos) + VOL4D*k];
//	   int z = get_spacecoord(spacepos, 3);
//	   correlator_ps[z] += spinor_squarenorm(tmp);
//	 }
//     }
//     }
//     // print the correlator in the "physical" normalization,
//    // i.e. in the mass basis, without kappa
//     // thus we have to multiply by (2*kappa)**2
//     printf("ps correlator:\n");
//     for(int i = 0; i<NSPACE; i++){
//       printf("%i\t(%.12e)\n", i, 4.*KAPPA*KAPPA*correlator_ps[i]);
//     }
//   }
//
//}


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

