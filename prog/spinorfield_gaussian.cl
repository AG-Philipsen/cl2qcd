__kernel void generate_gaussian_spinorfield(__global spinorfield * out, __global rngStateStorageType * rngStates)
{
	int local_size = get_local_size(0);
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	int loc_idx = get_local_id(0);
	int num_groups = get_num_groups(0);
	int group_id = get_group_id (0);
	int n, t;
	hmc_complex tmp;
	//sigma has to be 0.5 here
	hmc_float sigma = 0.5;
	spinor out_tmp;

#ifdef _SAME_RND_NUMBERS_
	if(id > 0) return;
	global_size = 1;
#endif

	hmc_ocl_ran rnd = loadRngState(rngStates);

	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size) {
		/** @todo this must be done more efficient */
		st_index pos = (id_tmp < VOLSPACE * NTIME / 2) ? get_even_site(id_tmp) : get_odd_site(id_tmp - (VOLSPACE * NTIME / 2));

		//CP: there are 12 complex elements in the spinor
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e0.e0.re = tmp.re;
		out_tmp.e0.e0.im = tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e0.e1.re = tmp.re;
		out_tmp.e0.e1.im = tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e0.e2.re = tmp.re;
		out_tmp.e0.e2.im = tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e1.e0.re = tmp.re;
		out_tmp.e1.e0.im = tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e1.e1.re = tmp.re;
		out_tmp.e1.e1.im = tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e1.e2.re = tmp.re;
		out_tmp.e1.e2.im = tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e2.e0.re = tmp.re;
		out_tmp.e2.e0.im = tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e2.e1.re = tmp.re;
		out_tmp.e2.e1.im = tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e2.e2.re = tmp.re;
		out_tmp.e2.e2.im = tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e3.e0.re = tmp.re;
		out_tmp.e3.e0.im = tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e3.e1.re = tmp.re;
		out_tmp.e3.e1.im = tmp.im;
		tmp = gaussianNormalPair(&rnd);
		out_tmp.e3.e2.re = tmp.re;
		out_tmp.e3.e2.im = tmp.im;

		//multiply by sigma
		out_tmp = real_multiply_spinor(out_tmp, sqrt(sigma));

		put_spinor_to_field(out_tmp, out, pos.space, pos.time);
	}

	storeRngState(rngStates, rnd);
}
