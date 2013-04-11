__kernel void set_zero_spinorfield_eoprec( __global spinorStorageType * const restrict x )
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);

	for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {
		putSpinor_eo(x, id_mem, set_spinor_zero());
	}
	return;
}
