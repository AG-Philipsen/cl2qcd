__kernel void create_volume_source(__global spinor * const restrict b, __global rngStateStorageType * const restrict rngStates)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

#ifdef _SAME_RND_NUMBERS_
	if(id > 0) return;
	global_size = 1;
#endif

	prng_state rnd;
	prng_loadState(&rnd, rngStates);

	/** @todo What is the correct norm here? */
	hmc_float sigma = 1. / ( VOL4D * 12. );
	spinor out_tmp;

	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size) {
		/** @todo this might be done more efficient */
		st_index pos = (id_tmp < VOLSPACE * NTIME / 2) ? get_even_site(id_tmp) : get_odd_site(id_tmp - (VOLSPACE * NTIME / 2));

		out_tmp = set_spinor_cold();

		//multiply by sigma
		out_tmp = real_multiply_spinor(out_tmp, sqrt(sigma));

		put_spinor_to_field(out_tmp, b, pos.space, pos.time);
	}

	prng_storeState(rngStates, &rnd);
}



