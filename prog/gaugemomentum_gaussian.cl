/** @todo rewrite gaussianComplexVector to handle structs??*/
__kernel void generate_gaussian_gaugemomenta(__global ae * out, __global hmc_ocl_ran * rnd)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);

	hmc_complex tmp;

#ifdef _SAME_RND_NUMBERS_
       if(id>0) return;
       global_size = 1;
#endif

	for(int id_tmp = id; id_tmp < GAUGEMOMENTASIZE; id_tmp += global_size) {
		//CP: there are 8 elements in ae
		tmp = gaussianNormalPair(&rnd[id]);
		out[id_tmp].e0 = tmp.re;
		out[id_tmp].e1 = tmp.im;
		tmp = gaussianNormalPair(&rnd[id]);
		out[id_tmp].e2 = tmp.re;
		out[id_tmp].e3 = tmp.im;
		tmp = gaussianNormalPair(&rnd[id]);
		out[id_tmp].e4 = tmp.re;
		out[id_tmp].e5 = tmp.im;
		tmp = gaussianNormalPair(&rnd[id]);
		out[id_tmp].e6 = tmp.re;
		out[id_tmp].e7 = tmp.im;
	}
}
