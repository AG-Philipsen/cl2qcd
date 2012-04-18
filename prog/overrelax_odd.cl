/** @file
 * Device code for the heatbath overrelaxation
 */

__kernel void overrelax_odd(__global Matrixsu3StorageType * gaugefield, const int mu, __global rngStateStorageType * rngStates)
{
	int id;

	hmc_ocl_ran rnd = loadRngState(rngStates);

	PARALLEL_FOR(id, VOLSPACE * NTIME / 2) {
		st_index pos = get_odd_site(id);
		perform_overrelaxing(gaugefield, mu, &rnd, pos.space, pos.time);
	}

	storeRngState(rngStates, rnd);
}
