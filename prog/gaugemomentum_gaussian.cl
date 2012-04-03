/** @todo rewrite gaussianComplexVector to handle structs??*/
__kernel void generate_gaussian_gaugemomenta(__global aeStorageType * const restrict out, __global rngStateStorageType * const restrict rngStates)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);

	hmc_ocl_ran rnd = loadRngState(rngStates);

	hmc_complex tmp;
#ifdef _SAME_RND_NUMBERS_
	if(id > 0) return;
	global_size = 1;
#endif
	for(int id_tmp = id; id_tmp < GAUGEMOMENTASIZE; id_tmp += global_size) {
		//CP: THERE ARE 8 ELEMENTS IN AE
		ae new_ae;

		tmp = gaussianNormalPair(&rnd);
		new_ae.e0 = tmp.re;
		new_ae.e1 =  tmp.im;
		tmp = gaussianNormalPair(&rnd);
		new_ae.e2 =  tmp.re;
		new_ae.e3 =  tmp.im;
		tmp = gaussianNormalPair(&rnd);
		new_ae.e4 =  tmp.re;
		new_ae.e5 =  tmp.im;
		tmp = gaussianNormalPair(&rnd);
		new_ae.e6 =  tmp.re;
		new_ae.e7 =  tmp.im;

		putAe(out, id_tmp, new_ae);
	}

	storeRngState(rngStates, rnd);
}
