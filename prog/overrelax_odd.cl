/** @file
 * Device code for the heatbath overrelaxation
 */

__kernel void overrelax_odd(__global ocl_s_gaugefield* gaugefield, const int mu, __global rngStateStorageType * rngStates)
{
	int id, id_tmp, size;
	id_tmp = get_global_id(0);
	size = get_global_size(0);

	hmc_ocl_ran rnd = loadRngState(rngStates);

	for(id = id_tmp; id < VOLSPACE * NTIME / 2; id += size) {
		st_index pos = get_odd_site(id);
		perform_overrelaxing(gaugefield, mu, &rnd, pos.space, pos.time, id_tmp);
	}

	storeRngState(rngStates, rnd);
}
