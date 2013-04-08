/** @todo rewrite gaussianComplexVector to handle structs??*/
__kernel void generate_gaussian_gaugemomenta(__global aeStorageType * const restrict out, __global rngStateStorageType * const restrict rngStates)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);

	prng_state rnd;
	prng_loadState(&rnd, rngStates);

	hmc_complex tmp;
#ifdef _SAME_RND_NUMBERS_
	if(id > 0) return;
	global_size = 1;
#endif
	PARALLEL_FOR(id_local, VOL4D_LOCAL) {
		site_idx site = get_site_idx(id_local >= (VOL4D_LOCAL / 2) ? get_even_st_idx_local(id_local - (VOL4D_LOCAL / 2)) : get_odd_st_idx_local(id_local));
		for(uint d = 0; d < 4; ++d) {
			const link_idx id_mem = get_link_idx(d, get_st_idx_from_site_idx(site));
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

			putAe(out, id_mem, new_ae);
		}
	}

	prng_storeState(rngStates, &rnd);
}
