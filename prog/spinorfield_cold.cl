/** @todo: here one can insert multi-threading, but this is not so important */
__kernel void set_spinorfield_cold(__global spinorfield* in)
{
	int id = get_global_id(0);
	if(id == 0) {
		hmc_float norm = 1. / sqrt((hmc_float)(12.f * VOLSPACE * NTIME));
		for(int t = 0; t < NTIME; t++) {
			for(int n = 0; n < VOLSPACE; n++) {
				in[get_global_pos(n, t)] = set_spinor_cold();
				in[get_global_pos(n, t)] = real_multiply_spinor(in[get_global_pos(n, t)], norm);
			}
		}
	}
}
