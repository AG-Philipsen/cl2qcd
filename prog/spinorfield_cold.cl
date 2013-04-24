/** @todo: here one can insert multi-threading, but this is not so important */
__kernel void set_spinorfield_cold(__global spinor * const restrict out)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	for(int id_mem = id; id_mem < SPINORFIELDSIZE_MEM; id_mem += global_size) {
		hmc_float norm = 1. / sqrt((hmc_float)(12.f * VOL4D_GLOBAL));
		out[id_mem] = real_multiply_spinor(set_spinor_cold(), norm);
	}

}
