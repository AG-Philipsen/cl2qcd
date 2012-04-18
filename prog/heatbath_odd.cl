__kernel void heatbath_odd(__global Matrixsu3StorageType * gaugefield, const int mu, __global rngStateStorageType * rngStates)
{
	int id, id_tmp, size;
	id_tmp = get_global_id(0);
	size = get_global_size(0);

	hmc_ocl_ran rnd = loadRngState(rngStates);

	for(id = id_tmp; id < VOLSPACE * NTIME / 2; id += size) {
		st_index pos = get_odd_site(id);
		perform_heatbath(gaugefield, mu, &rnd, pos.space, pos.time);
	}

	storeRngState(rngStates, rnd);
}
