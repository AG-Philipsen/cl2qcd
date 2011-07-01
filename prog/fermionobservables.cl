/**
 @file fermionmatrix
*/

//this is the pseudoscalar pion correlator in z-direction
__kernel void ps_correlator(__global spinorfield* phi){
   int local_size = get_local_size(0);
   int global_size = get_global_size(0);
   int id = get_global_id(0);
   int loc_idx = get_local_id(0);
   int num_groups = get_num_groups(0);
   int group_id = get_group_id (0);

   if(id==0){
	hmc_float correlator_ps[NSPACE];
	for(int i = 0; i<NSPACE; i++){
		correlator_ps[i]=0.;
	}
	for(int timepos = 0; timepos<NTIME; timepos++) {
   		for(int spacepos = 0; spacepos<VOLSPACE; spacepos++) {
		//correlator_ps[z] += |phi(n,t)|^2
		spinor tmp = phi[get_global_pos(spacepos, timepos)];
		int z = get_spacecoord(spacepos, 3);
		correlator_ps[z] += spinor_squarenorm(tmp);
   }}
	printf("ps correlator:\n");
	for(int i = 0; i<NSPACE; i++){
		printf("%i\t(%.12e)\n", i, correlator_ps[i]);
	}
  }

}