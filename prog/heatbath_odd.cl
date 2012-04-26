__kernel void heatbath_odd(__global Matrixsu3StorageType * const restrict dest_gaugefield, __global const Matrixsu3StorageType * const restrict src_gaugefield, const int mu, __global rngStateStorageType * const restrict rngStates)
{
	hmc_ocl_ran rnd;
	loadRngState(&rnd, rngStates);

	PARALLEL_FOR(id, VOLSPACE * NTIME / 2) {
		st_index pos = get_odd_site(id);
		perform_heatbath(dest_gaugefield, src_gaugefield, mu, &rnd, pos.space, pos.time);
	}

	storeRngState(rngStates, &rnd);
}

__kernel void heatbath_odd_hack(__global Matrixsu3StorageType * const restrict dest_gaugefield, __global const Matrixsu3StorageType * const restrict src_gaugefield, const int mu)
{
	PARALLEL_FOR(id, VOLSPACE * NTIME / 2) {
		st_index pos = get_odd_site(id);
		Matrixsu3 tmp = get_matrixsu3(src_gaugefield, pos.space, pos.time, mu);
		put_matrixsu3(dest_gaugefield, tmp, pos.space, pos.time, mu);
	}
}
