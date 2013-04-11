__kernel void heatbath_odd(__global Matrixsu3StorageType * const restrict gaugefield, const int mu, __global rngStateStorageType * const restrict rngStates)
{
	prng_state rnd;
	prng_loadState(&rnd, rngStates);

	PARALLEL_FOR(id, VOL4D_LOCAL / 2) {
		st_index pos = get_odd_st_idx_local(id);
		perform_heatbath(gaugefield, mu, &rnd, pos.space, pos.time);
	}

	prng_storeState(rngStates, &rnd);
}
