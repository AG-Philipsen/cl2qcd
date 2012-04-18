__kernel void heatbath_even(__global Matrixsu3StorageType * gaugefield, const int mu, __global rngStateStorageType * rngStates)
{
	int id;

	hmc_ocl_ran rnd = loadRngState(rngStates);

	PARALLEL_FOR(id, VOLSPACE * NTIME / 2) {
		st_index pos = get_even_site(id);
		perform_heatbath(gaugefield, mu, &rnd, pos.space, pos.time);
	}

	storeRngState(rngStates, rnd);
}

