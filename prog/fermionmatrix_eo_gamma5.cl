__kernel void gamma5_eoprec(__global spinorStorageType * const restrict inout)
{
	int global_size = get_global_size(0);
	int id = get_global_id(0);
	spinor out_tmp;

	for(int id_tmp = id; id_tmp < EOPREC_SPINORFIELDSIZE; id_tmp += global_size) {

		out_tmp = getSpinor_eo(inout, id_tmp);
		out_tmp = gamma5_local(out_tmp);
		putSpinor_eo(inout, id_tmp, out_tmp);
	}
}
