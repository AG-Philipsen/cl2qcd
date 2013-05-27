__kernel void set_zero_spinorfield_stagg_eoprec( __global staggeredStorageType * const restrict x )
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);

	for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {
		put_su3vec_to field_eo(x, id_mem, set_su3vec_zero());
	}
	return;
}
