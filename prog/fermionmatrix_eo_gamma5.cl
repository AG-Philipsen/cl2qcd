__kernel void gamma5_eo(__global spinorStorageType * const restrict inout)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	spinor out_tmp;

	for(int id_mem = id; id_mem < EOPREC_SPINORFIELDSIZE_MEM; id_mem += global_size) {

		out_tmp = getSpinor_eo(inout, id_mem);
		out_tmp = gamma5_local(out_tmp);
		putSpinor_eo(inout, id_mem, out_tmp);
	}
}
