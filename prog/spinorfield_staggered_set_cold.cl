/** @todo: here one can insert multi-threading, but this is not so important */
__kernel void set_cold_spinorfield_stagg(__global su3vec * const restrict out)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	for(int id_mem = id; id_mem < SPINORFIELDSIZE_MEM; id_mem += global_size) {
		hmc_float norm = 1. / sqrt((hmc_float)(3.f * VOL4D_GLOBAL));
		out[id_mem] = su3vec_times_real(set_su3vec_cold(), norm);
	}

}
