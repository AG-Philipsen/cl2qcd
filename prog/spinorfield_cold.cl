/** @todo: here one can insert multi-threading, but this is not so important */
__kernel void set_spinorfield_cold(__global spinorfield* in)
{
	int id = get_global_id(0);
	int global_size = get_global_size(0);

	for(int id_tmp = id; id_tmp < SPINORFIELDSIZE; id_tmp += global_size)
	{
	    st_index pos = (id_tmp % 2 == 0) ? get_even_site(id_tmp /   2) : get_odd_site(id_tmp / 2);

	    hmc_float norm = 1. / sqrt((hmc_float)(12.f * VOL4D));
	    in[get_global_pos(pos.space, pos.time)] = set_spinor_cold();
	    in[get_global_pos(pos.space, pos.time)] = real_multiply_spinor(in[get_global_pos(pos.space, pos.time)], norm);
	}

}
