__kernel void create_timeslice_source(__global spinor * const restrict b, __global rngStateStorageType * const restrict rngStates, int const timeslice)
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

	for(int id_tmp = id; id_tmp < VOLSPACE; id_tmp += global_size) {
		out_tmp = set_spinor_cold();

		//multiply by sigma
		out_tmp = real_multiply_spinor(out_tmp, sqrt(sigma));

		put_spinor_to_field(out_tmp, b, id_tmp, timeslice);
	}

	prng_storeState(rngStates, &rnd);

	return;
}
